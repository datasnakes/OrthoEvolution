# Standard Library
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# BioPython
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

from OrthoEvol.Orthologs.Align.guidance2 import Guidance2Commandline
from OrthoEvol.Orthologs.Align.orthoclustal import ClustalO
from OrthoEvol.Orthologs.Align.pal2nal import \
    PAL2NALCommandline as Pal2NalCommandline
from OrthoEvol.Orthologs.GenBank import GenBank
# OrthoEvol
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.utilities import FullUtilities


@dataclass(frozen=True, slots=True)
class _GuidancePaths:
    """Keep related GUIDANCE2 outputs together across filtering branches."""

    gene_directory: Path
    filtered_sequence: Path
    removed_sequences: Path
    alignment: Path
    sequence_column_filtered: Path
    column_filtered: Path
    masked: Path


class MultipleSequenceAlignment(object):
    """The MultipleSequenceAlignment (MSA) class uses the standard configuration
    along with function dispatching to give the end-user access to multiple
    alignment tools."""

    def __init__(self, project=None, project_path=os.getcwd(), genbank=GenBank, **kwargs):
        """Initialize the MultipleSequenceAlignment class.

        :param project: The project name.
        :type project: str or None
        :param project_path: The path to the project.
        :type project_path: str or Path
        :param genbank: The composer parameter which is used to configure the
                        GenBank class with the MSA class.
        :type genbank: class
        :param kwargs: The kwargs are used with the dispatcher as a way to
                        control the alignment pipeline. Can include:
                        - Guidance_config: Configuration for GUIDANCE2
                        - Pal2Nal_config: Configuration for PAL2NAL
                        - ClustalO_config: Configuration for ClustalO
        :type kwargs: dict
        :returns: If the kwargs are utilized with YAML or other
                  configurations, then this class returns an alignment
                  dictionary, which can be parsed to run specific alignment algorithms.
        """
        self.dispatcher_options = {"Guidance_config": ["GUIDANCE2", self.guidance2],
                                   "Pal2Nal_config": ["PAL2NAL", self.pal2nal],
                                   "ClustalO_config": ["CLUSTALO", self.clustalo]}
        # Set up loggers
        __log = LogIt()
        __logfile = None
        self.guidancelog = __log.default('guidance2', __logfile)
        self.pal2nallog = __log.default('pal2nal', __logfile)
        self.clustalolog = __log.default('clustalo', __logfile)

        # Initialize Utilities
        self.msa_utils = FullUtilities()

        # stop_codons = ['TAG', 'TAA', 'TGA']

        self.program = None
        self.alignment_dict = {}
        self.project = project
        self.project_path = project_path
        if project_path and project:
            self.project_path = Path(project_path) / Path(project)

        # Configuration of class attributes
        add_self = self.msa_utils.attribute_config(self, composer=genbank, checker=GenBank, project=project, project_path=project_path)
        for var, attr in add_self.__dict__.items():
            setattr(self, var, attr)

        # Determine which alignment to configure
        # And then run that alignment with the configuration.
        for config in self.dispatcher_options.keys():
            if config in kwargs.keys():
                program = self.dispatcher_options[config][0]
                aligner = self.dispatcher_options[config][1]
                aligner_configuration = kwargs[config]
                self.alignment_dict[program] = [aligner, aligner_configuration]

    def _build_guidance_paths(
        self,
        sequence_file: str | Path,
        sequence_type: str,
    ) -> _GuidancePaths:
        """Build every persistent output path before running GUIDANCE2."""
        suffixes = {
            "nuc": ("ffn", "na"),
            "aa": ("faa", "aa"),
        }
        try:
            sequence_extension, alignment_label = suffixes[sequence_type]
        except KeyError as error:
            supported_types = ", ".join(sorted(suffixes))
            raise ValueError(
                f"Unsupported GUIDANCE2 sequence type {sequence_type!r}; "
                f"expected one of: {supported_types}."
            ) from error

        gene = Path(sequence_file).stem
        gene_directory = Path(self.raw_data) / gene
        return _GuidancePaths(
            gene_directory=gene_directory,
            filtered_sequence=gene_directory / f"{gene}_G2.{sequence_extension}",
            removed_sequences=(
                gene_directory / f"{gene}_G2_removed.{sequence_extension}"
            ),
            alignment=gene_directory / f"{gene}_G2_{alignment_label}.aln",
            sequence_column_filtered=(
                gene_directory / f"{gene}_G2sfcf_{alignment_label}.aln"
            ),
            column_filtered=(
                gene_directory / f"{gene}_G2cf_{alignment_label}.aln"
            ),
            masked=gene_directory / f"{gene}_G2mf_{alignment_label}.aln",
        )

    def _run_guidance_command(self, **command_options: Any) -> None:
        """Run one configured command through the existing wrapper boundary."""
        command = Guidance2Commandline(**command_options)
        self.guidancelog.info(command)
        subprocess.check_call(
            [str(command)],
            stderr=subprocess.STDOUT,
            shell=True,
        )

    @staticmethod
    def _removed_sequence_file(iteration_directory: Path) -> Path:
        """Use the named GUIDANCE2 output when the program creates it."""
        named_file = (
            iteration_directory / "Seqs.Orig.fas.FIXED.Removed_Seq.With_Names"
        )
        if named_file.is_file():
            return named_file
        return iteration_directory / "Seqs.Orig.fas.FIXED.Removed_Seq"

    @staticmethod
    def _sequence_filter_output_directory(
        gene_directory: Path,
        column_filter: float | None,
        mask_filter: float | None,
    ) -> Path:
        """Select one stable directory for every filtering iteration."""
        if column_filter is not None:
            suffix = "sf_cf"
        elif mask_filter is not None:
            suffix = "sf_mf"
        else:
            suffix = "sf"
        return gene_directory / f"GUIDANCE2_{suffix}"

    @staticmethod
    def _should_stop_guidance(
        removed_sequence_count: int,
        iteration: int,
        maximum_iterations: int,
    ) -> bool:
        """Stop when GUIDANCE2 converges or reaches its configured limit."""
        return removed_sequence_count == 0 or iteration >= maximum_iterations

    def _run_guidance_sequence_filter(
        self,
        sequence_file: str | Path,
        msa_program: str,
        sequence_type: str,
        dataset: str,
        sequence_filter: str,
        column_filter: float | None,
        mask_filter: float | None,
        paths: _GuidancePaths,
        command_options: dict[str, Any],
        maximum_iterations: int | None,
        increment: float | None,
    ) -> Path:
        """Iteratively remove low-confidence sequences and return the last run."""
        if sequence_filter not in {"inclusive", "exclusive"}:
            raise ValueError(
                "sequence_filter must be either 'inclusive' or 'exclusive'."
            )

        if not isinstance(maximum_iterations, int) or maximum_iterations < 1:
            raise ValueError("iterations must be a positive integer.")
        if maximum_iterations > 1 and increment is None:
            raise ValueError("increment is required when iterations is greater than 1.")

        output_directory = self._sequence_filter_output_directory(
            paths.gene_directory,
            column_filter,
            mask_filter,
        )
        output_directory.mkdir(parents=True, exist_ok=True)

        for iteration in range(1, maximum_iterations + 1):
            iteration_directory = output_directory / f"iter_{iteration}"
            iteration_directory.mkdir(parents=True, exist_ok=True)

            if iteration > 1:
                if sequence_filter == "inclusive":
                    command_options["seqCutoff"] -= increment
                else:
                    command_options["seqCutoff"] += increment

            iteration_input = (
                Path(sequence_file) if iteration == 1 else paths.filtered_sequence
            )
            self._run_guidance_command(
                seqFile=str(iteration_input),
                msaProgram=msa_program,
                seqType=sequence_type,
                outDir=str(iteration_directory),
                **command_options,
            )

            removed_file = self._removed_sequence_file(iteration_directory)
            removed_sequence_count = sum(
                1 for _ in SeqIO.parse(str(removed_file), "fasta")
            )

            if iteration == 1:
                SeqIO.write(
                    SeqIO.parse(str(removed_file), "fasta"),
                    str(paths.removed_sequences),
                    "fasta",
                )
            elif removed_sequence_count > 0:
                self.msa_utils.multi_fasta_manipulator(
                    str(paths.removed_sequences),
                    str(removed_file),
                    str(paths.removed_sequences),
                    manipulation="add",
                )

            if removed_sequence_count > 0:
                self.msa_utils.multi_fasta_manipulator(
                    str(sequence_file),
                    str(paths.removed_sequences),
                    str(paths.filtered_sequence),
                    manipulation="remove",
                )

            if self._should_stop_guidance(
                removed_sequence_count,
                iteration,
                maximum_iterations,
            ):
                filtered_alignment = (
                    iteration_directory
                    / f"{dataset}.{msa_program}.aln.Sorted.With_Names"
                )
                renamed_alignment = shutil.copy(
                    str(filtered_alignment),
                    str(paths.alignment),
                )
                self.msa_utils.multi_fasta_manipulator(
                    str(renamed_alignment),
                    str(sequence_file),
                    str(renamed_alignment),
                    manipulation="sort",
                )
                return iteration_directory

        raise RuntimeError("GUIDANCE2 sequence filtering ended without a result.")

    def _apply_guidance_post_filter(
        self,
        sequence_file: str | Path,
        msa_program: str,
        sequence_type: str,
        dataset: str,
        column_filter: float | None,
        mask_filter: float | None,
        paths: _GuidancePaths,
        iteration_directory: Path,
    ) -> None:
        """Apply the optional column or residue filter after sequence filtering."""
        if column_filter is not None:
            filtered_alignment = (
                iteration_directory
                / f"{dataset}.{msa_program}.Without_low_SP_Col.With_Names"
            )
            shutil.copy(
                str(filtered_alignment),
                str(paths.sequence_column_filtered),
            )
        elif mask_filter is not None:
            alignment_to_mask = (
                iteration_directory / f"{dataset}.{msa_program}.aln.With_Names"
            )
            residue_pair_scores = (
                iteration_directory
                / f"{dataset}.{msa_program}.Guidance2_res_pair_res.scr"
            )
            self._run_guidance_command(
                align=False,
                seqType=sequence_type,
                maskCutoff=mask_filter,
                maskFile=str(alignment_to_mask),
                rprScores=str(residue_pair_scores),
                output=str(paths.masked),
            )
            self.msa_utils.multi_fasta_manipulator(
                str(paths.masked),
                str(sequence_file),
                str(paths.masked),
                manipulation="sort",
            )

    def guidance2(
        self,
        seqFile: str | Path,
        msaProgram: str,
        seqType: str,
        dataset: str = "MSA",
        seqFilter: str | None = None,
        columnFilter: float | None = None,
        maskFilter: float | None = None,
        **kwargs: Any,
    ) -> None:
        """Run the GUIDANCE2 command line wrapper from BioPython.

        The Guidance2 algorithm is used to filter sequence alignments in
        different ways.  Here we employ a few of our own strategies on top of
        Guidance2.

        :param seqFile:  The sequence file required by GUIDANCE2.
        :param msaProgram:  The msa program to be used by GUIDANCE2.
                            ("CLUSTALW", "PRANK", "MAFFT", or "MUSCLE")
        :param seqType: The type of sequences to be aligned in GUIDANCE2.
                        ("aa", "nuc", or "codon")
        :param dataset: The name of the dataset, which is used for file
                        naming convention among other things in GUIDANCE2.
        :param seqFilter: The sequence filter parameter is None, "inclusive",
                          or "exclusive".  If inclusive the SeqCutoff
                          decreases for every iteration.  If exclusive the
                          SeqCutoff increases for every iteration, and so the
                          algorithm excludes more genes from the alignment.
                          (An OrthoEvol strategy)
        :param columnFilter: The column filter removes columns from the
                             alignment using GUIDANCE2.
        :param maskFilter: The mask filter uses GUIDANCE2 maskLowScoresResidue
                           script to mask the low scoring residues.
        :param kwargs: The kwargs are used to configure GUIDANCE2 with
                       specific parameters including seqCutoff and colCutoff.
                       It can also be used to set the number of iterations and
                       the increment number, which controls how seqCutoff and
                       colCutoff change for each iteration.
        :return:  Returns Guidance2 files.
        """

        if columnFilter is not None and maskFilter is not None:
            raise ValueError("columnFilter and maskFilter are mutually exclusive.")

        self.guidancelog.info("Guidance2 will be used.")
        self.program = "GUIDANCE2"
        paths = self._build_guidance_paths(seqFile, seqType)
        self.guidancelog.info(paths.gene_directory)

        command_options = dict(kwargs)
        maximum_iterations = command_options.pop("iterations", None)
        increment = command_options.pop("increment", None)
        command_options.setdefault("seqCutoff", 0.6)
        command_options.setdefault("colCutoff", 0.93)

        if seqFilter is not None:
            iteration_directory = self._run_guidance_sequence_filter(
                seqFile,
                msaProgram,
                seqType,
                dataset,
                seqFilter,
                columnFilter,
                maskFilter,
                paths,
                command_options,
                maximum_iterations,
                increment,
            )
            self._apply_guidance_post_filter(
                seqFile,
                msaProgram,
                seqType,
                dataset,
                columnFilter,
                maskFilter,
                paths,
                iteration_directory,
            )
        elif columnFilter is not None:
            output_directory = paths.gene_directory / "GUIDANCE2_cf"
            output_directory.mkdir(parents=True, exist_ok=True)
            self._run_guidance_command(
                seqFile=str(seqFile),
                msaProgram=msaProgram,
                seqType=seqType,
                outDir=str(output_directory),
                **command_options,
            )
            filtered_alignment = (
                output_directory
                / f"{dataset}.{msaProgram}.Without_low_SP_Col.With_Names"
            )
            shutil.copy(str(filtered_alignment), str(paths.column_filtered))
        elif maskFilter is not None:
            required_options = {"aln2mask", "rprScores", "maskedFile"}
            missing_options = required_options.difference(command_options)
            if missing_options:
                missing = ", ".join(sorted(missing_options))
                raise ValueError(f"Missing mask options: {missing}.")

            masked_file = command_options.pop("maskedFile")
            self._run_guidance_command(
                align=False,
                seqType=seqType,
                maskCutoff=maskFilter,
                maskFile=command_options.pop("aln2mask"),
                rprScores=command_options.pop("rprScores"),
                output=masked_file,
            )
            self.msa_utils.multi_fasta_manipulator(
                masked_file,
                str(seqFile),
                masked_file,
                manipulation="sort",
            )

    def pal2nal(self, aa_alignment, na_fasta, output_type='paml', nogap=True, nomismatch=True, downstream='paml'):
        """This Pal2Nal method works with the Pal2Nal command line wrapper.

        It uses a protein alignment to generate a codon alignment from the
        corresponding nucleic acid sequences.  This is useful for downstream
        PAML analysis. This function also catches and removes taxa that are
        inconsistent with Pal2Nal's algorithm.

        :param aa_alignment: An amino acid alignment that is used as a guide
                             for a nucleic acid alignment.
        :param na_fasta: The FASTA file that contains matching/ordered
                         sequences corresponding to the aa_alignment.
        :param output_type: The format of the resulting alignment.
                            ("clustal", "paml", "fasta", "codon")
        :param nogap: Removes the gaps and in-frame stop codons from the
                      alignment to work better with PAML.
        :param nomismatch: Removes mismatched codons between protein and DNA
                           sequences.
        :param downstream: Used as a naming convention for a better and more
                           obvious pipeline.
        :return: A codon alignment.
        """

        removed = []
        # Create output directory for PAL2NAL
        outDir = 'PAL2NAL'
        gene = Path(na_fasta).stem
        geneDir = self.raw_data / Path(gene)
        outDir = geneDir / Path(outDir)
        Path.mkdir(outDir, exist_ok=True)

        # Name the output file using the gene and the downstream application
        output_file = str(geneDir / Path(gene + '_P2N.%s.aln' % downstream))

        # Create an alignment
        P2Ncmd = Pal2NalCommandline(pepaln=aa_alignment, nucfasta=na_fasta, output_file=output_file, output=output_type,
                                    nogap=nogap, nomismatch=nomismatch)
        self.pal2nallog.info(P2Ncmd)

        # Use a while loop to catch errors and remove sequences that aren't working with pal2nal
        pal2nal_flag = True
        while pal2nal_flag is True:
            pal2nal = subprocess.Popen([str(P2Ncmd)], stderr=subprocess.PIPE,
                                       stdout=subprocess.PIPE, shell=True,
                                       encoding='utf-8')
            error = pal2nal.stderr.readlines()
            out = pal2nal.stdout.readlines()
            pal2nal.wait()

            # Catch errors
            if 'ERROR: inconsistency between the following pep and nuc seqs' in error[0]:
                self.pal2nallog.warning('Caught the pal2nal error!')
                self.pal2nallog.warning(error[0])
                for err in error:
                    if '>' in err:
                        removed.append(err.strip('>' '\n'))
                self.msa_utils.multi_fasta_manipulator(na_fasta, removed, na_fasta)
                self.msa_utils.multi_fasta_manipulator(aa_alignment, removed, aa_alignment)

            # If no errors then break the while loop
            else:
                if removed is not None:
                    # If any sequences were removed then write them to a file.
                    p2n_remFile = str(geneDir / Path(gene + '_P2N_removed.txt'))
                    with open(p2n_remFile, 'w') as p2n_rem:
                        for name in removed:
                            p2n_rem.write(name)
                pal2nal_flag = False

            self.pal2nallog.info('Error: ' + str(error))
            self.pal2nallog.info('Out: ' + str(out))

    def clustalo(self, infile, outfile, outfmt="fasta"):
        """Align protein/amino acid sequences using Clustal Omega.

        :param infile: Input a multifasta protein file.
        :param outfile: Output an aligned multifasta file.
        :param outfmt:  (Default value = "fasta")
        """
        try:
            clustalo = ClustalO(infile=infile, outfile=outfile, outfmt=outfmt)
            clustalo.runclustalomega()
        except Exception:
            pass
