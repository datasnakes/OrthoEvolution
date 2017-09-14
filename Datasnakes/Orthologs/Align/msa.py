import os
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

from Datasnakes.Orthologs.Align.guidance2 import Guidance2Commandline
from Datasnakes.Orthologs.Align.pal2nal import Pal2NalCommandline
from Datasnakes.Orthologs.Genbank.genbank import GenBank
from Datasnakes.Orthologs.Genbank.utils import multi_fasta_manipulator


class MultipleSequenceAlignment(object):

    def __init__(self, project, aln_program, project_path=None, genbank=GenBank, **kwargs):

        self.program = aln_program
        self.project = project
        if isinstance(genbank, GenBank):
            setattr(genbank, 'project', project)
            for key, value in genbank.__dict__.items():
                setattr(self, key, value)
            print('project_path=%s' % self.project_path)
        else:
            if project_path:
                self.project_path = Path(project_path) / Path(self.project)
            else:
                self.project_path = Path(os.getcwd()) / Path(self.project)
            Path.mkdir(self.project_path, parents=True, exist_ok=True)
            print('project_path=%s' % self.project_path)
            self.removed_gb_config(kwargs)

        if kwargs:
            print('aln-kwargs')
            if aln_program == 'GUIDANCE2':
                self.align = self.guidance2
                self.guidance2(**kwargs)
            elif aln_program == 'CLUSTALO':
                self.align = self.clustalo
                self.clustalo(**kwargs)
            elif aln_program == 'PAL2NAL':
                self.align = self.pal2nal
                self.pal2nal(**kwargs)
        else:
            print('aln-notkwargs')
            if aln_program is 'GUIDANCE2':
                self.align = self.guidance2
                print(self.align)
            elif aln_program is 'CLUSTALO':
                self.align = self.clustalo
                print(self.align)
            elif aln_program is 'PAL2NAL':
                self.align = self.pal2nal
                print(self.align)


    def removed_gb_config(self, kwargs):
        self.raw_data = self.project_path / Path('raw_data')
        self.data = self.project_path / Path('data')

        for key, value in kwargs.items():
            setattr(self, key, value)

        Path.mkdir(self.raw_data, exist_ok=True)
        Path.mkdir(self.data, exist_ok=True)

    def guidance2(self, seqFile, msaProgram, seqType, dataset='MSA', seqFilter=None, columnFilter=None, maskFilter=None, **kwargs):
        # Name and Create the output directory
        outDir = self.program
        gene = Path(seqFile).stem
        geneDir = self.raw_data / Path(gene)

        if seqType is 'nuc':
            g2_seqFile = str(geneDir / Path(gene + '_G2.ffn'))  # Need for all iterations
            rem_file = str(geneDir / Path(gene + '_G2_removed.ffn'))   # Need for all iterations
            g2_alnFile = str(geneDir / Path(gene + '_G2_na.aln'))
            g2_seqcolFilter = str(geneDir / Path(gene + 'G2sfcf_na.aln'))
            g2_colFilter = str(geneDir / Path(gene + '_G2cf_na.aln'))
            g2_maskedFile = str(geneDir / Path(gene + '_G2mf_na.aln'))
        elif seqType is 'aa':
            g2_seqFile = str(geneDir / Path(gene + '_G2.faa'))  # Need for all iterations
            rem_file = str(geneDir / Path(gene + '_G2_removed.faa'))   # Need for all iterations
            g2_alnFile = str(geneDir / Path(gene + '_G2_aa.aln'))
            g2_colFilter = str(geneDir / Path(gene + '_G2cf_aa.aln'))
            g2_maskedFile = str(geneDir / Path(gene + '_G2mf_aa.aln'))

        # Add the Guidance 2 cutoffs to the keyword arguments
        if 'seqCutoff' not in kwargs.keys():
            kwargs['seqCutoff'] = 0.6
        if 'colCutoff' not in kwargs.keys():
            kwargs['colCutoff'] = 0.93

        # Filter Sequences, and then either remove columns, mask residues or do nothing
        if seqFilter is not None:

            # Filter "bad" sequences and iterate over the good ones until the cutoff is reached.
            iterFlag = True
            iteration = 0
            while iterFlag is True:
                set_iter = kwargs['iterations']
                iteration += 1
                # Create paths for output files
                if columnFilter is not None:
                    outDir = self.raw_data / Path(gene) / Path(outDir + '_sf_cf')
                elif maskFilter is not None:
                    outDir = self.raw_data / Path(gene) / Path(outDir + '_sf_mf')
                else:
                    outDir = self.raw_data / Path(gene) / Path(outDir + '_sf')
                Path.mkdir(outDir, parents=True, exist_ok=True)
                iterDir = Path(outDir) / Path('iter_%s' % iteration)
                g2_rem_file = str(iterDir / Path('Seqs.Orig.fas.FIXED.Removed_Seq.With_Names'))  # Need for all iterations
                Path.mkdir(iterDir, parents=True, exist_ok=True)

                # Create files for masking
                if maskFilter is not None:
                    g2_aln2mask = str(iterDir / Path('%s.%s.aln.With_Names' % (dataset, msaProgram)))
                    g2_rprScores = str(iterDir / Path('%s.%s.Guidance2_res_pair_res.scr' % (dataset, msaProgram)))

                if iteration == 1:

                    # seqFile is the given input
                    G2Cmd = Guidance2Commandline(seqFile=seqFile, msaProgram=msaProgram, seqType=seqType,
                                                 outDir=str(iterDir), **kwargs)
                    print(G2Cmd)
                    subprocess.check_call([str(G2Cmd)], stderr=subprocess.STDOUT, shell=True)
                    # Copy the Guidance removed seq file and paste it to the home directory
                    # Creates the rem_file
                    # Files without any removed don't have the file *.With_Names
                    if os.path.isfile(g2_rem_file) is False:
                        g2_rem_file = str(iterDir / Path('Seqs.Orig.fas.FIXED.Removed_Seq'))
                    SeqIO.write(SeqIO.parse(g2_rem_file, 'fasta'), rem_file, 'fasta')  # Need for iter_1

                    # Filter the input NA fasta file using Guidance output
                    # Creates the g2_seqFile
                    multi_fasta_manipulator(seqFile, g2_rem_file, g2_seqFile, manipulation='remove')  # Do after copying (iter_1) or adding (iter_n)
                    iterFlag = True

                elif set_iter >= iteration > 1:

                    # Depending on the filter strategy increment the seqCutoff
                    if seqFilter is "inclusive":
                        kwargs['seqCutoff'] -= kwargs['increment']
                    elif seqFilter is "exclusive":
                        kwargs['seqCutoff'] += kwargs['increment']
                    # seqFile changes to g2_seqFile and the cutoffs change
                    G2Cmd = Guidance2Commandline(seqFile=g2_seqFile, msaProgram=msaProgram, seqType=seqType,
                                                 outDir=str(iterDir), **kwargs)
                    print(G2Cmd)
                    subprocess.check_call([str(G2Cmd)], stderr=subprocess.STDOUT, shell=True)

                    # Get the removed sequence count
                    rem_count = 0
                    # Files without any removed don't have the file *.With_Names
                    if os.path.isfile(g2_rem_file) is False:
                        g2_rem_file = str(Path(iterDir) / Path('Seqs.Orig.fas.FIXED.Removed_Seq'))

                    for rec in SeqIO.parse(g2_rem_file, 'fasta'):  # Need for all iterations
                        rem_count += 1

                    # If sequences are removed, then iterate again on the "good" sequences
                    if rem_count > 0:
                        # Add new sequences to the rem_file
                        multi_fasta_manipulator(rem_file, g2_rem_file, rem_file, manipulation='add')
                        # Filter the input fasta file using the updated rem_file
                        multi_fasta_manipulator(seqFile, rem_file, g2_seqFile, manipulation='remove')
                        iterFlag = True
                    # If sequences aren't removed, then stop iterating
                    if rem_count < 0 or set_iter == iteration:
                        filtered_alignment = Path(iterDir) / Path('%s.%s.aln.Sorted.With_Names' % (dataset, msaProgram))
                        renamed_alignment = shutil.copy(str(filtered_alignment), g2_alnFile)
                        multi_fasta_manipulator(str(renamed_alignment), str(seqFile), str(renamed_alignment), manipulation='sort')
                        iterFlag = False

            if columnFilter is not None:
                col_filt_align = iterDir / Path('%s.%s.Without_low_SP_Col.With_Names' % (dataset, msaProgram))
                shutil.copy(str(col_filt_align), g2_seqcolFilter)

            elif maskFilter is not None:
                G2Cmd = Guidance2Commandline(align=False, seqFile=seqFile, msaProgram=msaProgram, seqType=seqType,
                                             outDir=str(iterDir), maskCutoff=maskFilter, maskFile=g2_aln2mask,
                                             rprScores=g2_rprScores, output=g2_maskedFile, **kwargs)
                print(G2Cmd)
                subprocess.check_call([str(G2Cmd)], stderr=subprocess.STDOUT, shell=True)
                multi_fasta_manipulator(g2_maskedFile, str(seqFile), g2_maskedFile, manipulation='sort')

        # Only COLUMN FILTER the bad columns
        elif columnFilter is not None:
            outDir = self.raw_data / Path(gene) / Path(outDir + '_cf')
            Path.mkdir(outDir, parents=True, exist_ok=True)
            G2Cmd = Guidance2Commandline(seqFile=seqFile, msaProgram=msaProgram, seqType=seqType,
                                         outDir=str(outDir), **kwargs)
            print(G2Cmd)
            subprocess.check_call([str(G2Cmd)], stderr=subprocess.STDOUT, shell=True)
            col_filt_align = outDir / Path('%s.%s.Without_low_SP_Col.With_Names' % (dataset, msaProgram))
            shutil.copy(str(col_filt_align), g2_colFilter)
            print()
        # Only MASK the bad residues
        elif maskFilter is not None:
            outDir = self.raw_data / Path(gene) / Path(outDir + '_sf')
            G2Cmd = Guidance2Commandline(seqFile=seqFile, msaProgram=msaProgram, seqType=seqType,
                                         outDir=str(outDir), maskCutoff=maskFilter, maskFile=kwargs['aln2mask'],
                                         rprScores=kwargs['rprScores'], output=kwargs['maskedFile'], **kwargs)
            print(G2Cmd)
            subprocess.check_call([str(G2Cmd)], stderr=subprocess.STDOUT, shell=True)
            multi_fasta_manipulator(kwargs['maskedFile'], str(seqFile), kwargs['maskedFile'], manipulation='sort')

    def pal2nal(self, aa_alignment, na_fasta, output_type='paml', nogap=True, nomismatch=True, downstream='paml'):
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
        print(P2Ncmd)

        # Use a while loop to catch errors and remove sequences that aren't working with pal2nal
        pal2nal_flag = True
        while pal2nal_flag is True:
            pal2nal = subprocess.Popen([str(P2Ncmd)], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True,
                                       encoding='utf-8')
            error = pal2nal.stderr.readlines()
            out = pal2nal.stdout.readlines()
            pal2nal.wait()

            # Catch errors
            if 'ERROR: inconsistency between the following pep and nuc seqs' in error[0]:
                print('Caught the pal2nal error!')
                print(error[0])
                for err in error:
                    if '>' in err:
                        removed.append(err.strip('>' '\n'))
                multi_fasta_manipulator(na_fasta, removed, na_fasta)
                multi_fasta_manipulator(aa_alignment, removed, aa_alignment)

            # If no errors then break the while loop
            else:
                if removed is not None:
                    # If any sequences were removed then write them to a file.
                    p2n_remFile = str(geneDir / Path(gene + '_P2N_removed.txt'))
                    with open(p2n_remFile, 'w') as p2n_rem:
                        for name in removed:
                            p2n_rem.write(name)
                pal2nal_flag = False

            print('Error: ' + str(error))
            print('Out: ' + str(out))

    def clustalo(self, infile, outfile, logpath, outfmt="fasta"):
        """This class aligns amino acids sequences using parameters similar to
        the default parameters.

        These parameters include 2 additional iterations for the hmm.
        """
        clustalo_cline = ClustalOmegaCommandline(infile=infile, cmd="clustalo",
                                                 outfile=outfile, seqtype="PROTEIN",
                                                 max_hmm_iterations=2, infmt="fasta",
                                                 outfmt=outfmt, iterations=3,
                                                 verbose=True,
                                                 force=True, log=logpath)
        stdout, stderr = clustalo_cline()

        # Run the command
        clustalo_cline()
        if stderr:
            print(stderr)
        if stdout:
            print(stdout)





