import os
from pathlib import Path
from shutil import copy
from Datasnakes.Orthologs.Align.QualityControl.guidance2 import Guidance2Commandline
from Datasnakes.Orthologs.Align.QualityControl.pal2nal import Pal2NalCommandline
from Bio import SeqIO
from Datasnakes.Orthologs.Genbank.utils import multi_fasta_manipulator


# For one gene at a time
class QCAlignmentCommandline(Guidance2Commandline, Pal2NalCommandline):

    def __init__(self, na_fasta, aa_fasta, home=os.getcwd(), msaProgram='CLUSTALW', na_bootstraps=10, aa_bootstraps=25,
                 na_seqCutoff=0.6, aa_seqCutoff=0.6, na_colCutoff=0, aa_colCutoff=0.88):
        Guidance2Commandline.__init__(outOrder='as_input', dataset=Path(na_fasta).stem, msaProgram=msaProgram)
        Pal2NalCommandline.__init__(nogap=True, nomismatch=True, output_file=Path(na_fasta).stem + '_P2N.aln')

        # Initialize
        self.home = Path(home)
        self.na_guidance_path = self.home / Path('NA_Guidance2')
        self.aa_guidance_path = self.home / Path('AA_Guidance2')
        self.gene = Path(self.na_fasta).stem

        # Guidance2 iterations over the nucleic acid sequences
        iteration_flag = True
        iteration = 1
        seqFile = Path(self.home) / Path(na_fasta)
        while iteration_flag is True:
            if iteration > 1:
                na_seqCutoff = 0.7
                na_bootstraps = 15
            iteration_flag, seqFile = self.nucleic_acid_guidance(iteration, seqFile=seqFile,
                                                                 outDir=self.na_guidance_path, bootstraps=na_bootstraps,
                                                                 seqCutoff=na_seqCutoff, colCutoff=na_colCutoff)
            iteration += 1
        # Guidance 2 amino acid alignment filter
        aa_alignment = self.amino_acid_guidance(seqFile=Path(self.home) / Path(aa_fasta), remFile=seqFile,
                                 outDir=self.aa_guidance_path, bootstraps=aa_bootstraps,
                                 seqCutoff=aa_seqCutoff, colCutoff=aa_colCutoff)
        na_fasta = self.home / Path(self.gene + '_G2.ffn')
        na_alignment = self.home / Path(self.gene + '_P2N_na.aln')
        # PAL2NAL nucleic acid alignment
        self.pal2nal_conversion(aa_alignment, na_fasta, na_alignment, )

    def nucleic_acid_guidance(self, iteration, seqFile, outDir, bootstraps, seqCutoff, colCutoff):
        seqType = 'nuc'

        outDir = Path(outDir) / 'iter_%s' % iteration
        Path.mkdir(outDir, parents=True, exist_ok=True)

        G2Cmd = Guidance2Commandline(seqFile=seqFile, seqType=seqType, outDir=outDir, bootstraps=bootstraps,
                                     seqCutoff=seqCutoff, colCutoff=colCutoff)

        G2Cmd()
        removed_file = Path(outDir) / Path('Seqs.Orig.fas.FIXED.Removed_Seq')
        renamed_file1 = copy(removed_file.__str__(), Path(self.gene + '_G2_removed.ffn').__str__())
        moved_file1 = copy(renamed_file1, self.home.__str__())

        seq_file = Path(outDir) / Path('Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names')
        renamed_file2 = copy(seq_file.__str__(), Path(self.gene + '_G2.ffn').__str__())
        moved_file2 = copy(renamed_file2.__str__(), self.home.__str__())

        # Copy and rename files
        rem_count = SeqIO.write(SeqIO.parse(removed_file, 'fasta'), renamed_file1, 'fasta')
        if iteration == 1:
            iter_flag = True
            return_file = moved_file2
        elif iteration > 1:
            if rem_count > 0:
                multi_fasta_manipulator(moved_file1, removed_file, manipulation='add', added_name='')
                iter_flag = True
                return_file = moved_file2
            else:
                iter_flag = False
                return_file = moved_file1
        return iter_flag, return_file

    def amino_acid_guidance(self, seqFile, remFile, outDir, bootstraps, seqCutoff, colCutoff):
        seqType = 'aa'
        filtered_fasta = multi_fasta_manipulator(seqFile, remFile)
        G2Cmd = Guidance2Commandline(seqFile=filtered_fasta, seqType=seqType, outDir=outDir, bootstraps=bootstraps,
                                     seqCutoff=seqCutoff, colCutoff=colCutoff)
        G2Cmd()

        filtered_alignment = Path(outDir) / Path('MSA.CLUSTALW.Without_low_SP_Col.With_Names')
        renamed_alignment = copy(filtered_alignment.__str__(), Path(self.gene + '_G2_aa.aln').__str__())
        moved_alignment = copy(renamed_alignment.__str__(), self.home.__str__())
        print('Align the filtered amino acid sequences using guidance 2')
        return Path(moved_alignment)

    def pal2nal_conversion(self, aa_alignment, na_fasta, output_file):

        P2Ncmd = Pal2NalCommandline(pepaln=aa_alignment, nucfasta=na_fasta, output_file=output_file, nogap=True,
                                    nomismatch=True)
        P2Ncmd()

        print('Align the nucleic acids using the amino acid alignment.')
