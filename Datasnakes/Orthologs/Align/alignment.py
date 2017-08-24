from Datasnakes.Orthologs.Align.QualityControl.guidance2 import Guidance2Commandline
from Datasnakes.Orthologs.Align.QualityControl.pal2nal import Pal2NalCommandline
from Datasnakes.Orthologs.Genbank.genbank import GenBank
from Datasnakes.Orthologs.Genbank.utils import multi_fasta_manipulator
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from pathlib import Path
import os
import subprocess
import shutil


class Alignment(GenBank):

    def __init__(self, program, **kwargs):
        super().__init__(**kwargs)

        self.program = program
        if kwargs:
            if program == 'GUIDANCE2':
                self.align = self.guidance2
                self.guidance2(kwargs)
            elif program == 'CLUSTALO':
                self.align = self.clustalo
                self.clustalo(kwargs)
            elif program == 'PAL2NAL':
                self.align = self.pal2nal
                self.pal2nal(kwargs)
        else:
            if program == 'GUIDANCE2':
                self.align = self.guidance2
            elif program == 'CLUSTALO':
                self.align = self.clustalo
            elif program == 'PAL2NAL':
                self.align = self.pal2nal

        print()

    def guidance2(self, seqFile, msaProgram, seqType, dataset='MSA', seqFilter=None, columnFilter=None, maskFilter=None, **kwargs):
        # Name and Create the output directory
        outDir = self.program
        gene = Path(seqFile).stem
        geneDir = self.raw_data / Path(gene)

        if seqType is 'nuc':
            g2_seqFile = str(geneDir / Path(self.gene + '_G2.ffn'))  # Need for all iterations
            rem_file = str(geneDir / Path(self.gene + '_G2_removed.ffn'))   # Need for all iterations
            g2_alnFile = str(geneDir / Path(self.gene + '_G2_na.aln'))
            g2_seqcolFilter = str(geneDir / Path(self.gene + 'G2sfcf_na.aln'))
            g2_colFilter = str(geneDir / Path(self.gene + '_G2cf_na.aln'))
            g2_maskFilter = str(geneDir / Path(self.gene + '_G2mf_na.aln'))
        elif seqType is 'aa':
            g2_seqFile = str(geneDir / Path(self.gene + '_G2.faa'))  # Need for all iterations
            rem_file = str(geneDir / Path(self.gene + '_G2_removed.faa'))   # Need for all iterations
            g2_alnFile = str(geneDir / Path(self.gene + '_G2_aa.aln'))
            g2_colFilter = str(geneDir / Path(self.gene + '_G2cf_aa.aln'))
            g2_maskFilter = str(geneDir / Path(self.gene + '_G2mf_aa.aln'))

        # Add the Guidance 2 cutoffs to the keyword arguments
        if 'seqCutoff' not in kwargs.keys():
            kwargs['seqCutoff'] = 0.6
        if 'colCutoff' not in kwargs.keys():
            kwargs['colCutoff'] = 0.93

        # filter = "inclusive", "exclusive", "neutral", or None (default) or "mask" (columnFilter)
        if seqFilter is not None:

            # Filter "bad" sequences and iterate over the good ones until the cutoff is reached.
            iterFlag = True
            iteration = 0
            while iterFlag is True:

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

                elif iteration > 1:

                    # Depending on the filter strategy increment the seqCutoff
                    if seqFilter is "inclusive":
                        kwargs['seqCutoff'] -= kwargs['increment']
                    elif seqFilter is "exclusive":
                        kwargs['seqCutoff'] += kwargs['increment']
                    # seqFile changes to g2_seqFile and the cutoffs change
                    G2Cmd = Guidance2Commandline(seqFile=g2_seqFile, msaProgram=msaProgram, seqType=seqType, outDir=str(iterDir), **kwargs)
                    print(G2Cmd)
                    subprocess.check_call([str(G2Cmd)], stderr=subprocess.STDOUT, shell=True)

                    # Get the removed sequence count
                    rem_count = 0
                    # Files without any removed don't have the file *.With_Names
                    if os.path.isfile(g2_rem_file) is False:
                        g2_rem_file = str(Path(iterDir) / Path('Seqs.Orig.fas.FIXED.Removed_Seq'))

                    for rec in SeqIO.parse(g2_rem_file, 'fasta'):  # Need for all iterations
                        rem_count += 1

                    if rem_count > 0:
                        # Add new sequences to the rem_file
                        multi_fasta_manipulator(rem_file, g2_rem_file, rem_file, manipulation='add')
                        # Filter the input fasta file using the updated rem_file
                        multi_fasta_manipulator(seqFile, rem_file, g2_seqFile, manipulation='remove')
                        iterFlag = True
                    else:
                        filtered_alignment = Path(iterDir) / Path('%s.%s.aln.Sorted.With_Names' % (dataset, msaProgram))
                        renamed_alignment = shutil.copy(str(filtered_alignment), g2_alnFile)
                        multi_fasta_manipulator(str(renamed_alignment), str(seqFile), str(renamed_alignment), manipulation='sort')
                        iterFlag = False
            if columnFilter is not None:
                col_filt_align = iterDir / Path('%s.%s.Without_low_SP_Col.With_Names' % (dataset, msaProgram))
                shutil.copy(str(col_filt_align), g2_seqcolFilter)
            elif maskFilter is not None:
                print()
                # TODO-ROB: create masking filter command line option

        # Filter the "bad" columns using good sequences
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
        # Mask the "bad" residues using the good sequences
        elif maskFilter is not None:
            print()

    def pal2nal(self, aa_alignment, na_fasta, output_file, nogap=True, nomismatch=True):
        # TODO-ROB:  Add a filter step so if error[0] = Error: #---  ERROR: inconsistency between the following pep and nuc seqs  ---#
        # TODO-ROB:  ....then remove error[2] which is the sequence it can't read
        removed = []
        # Create output directory for PAL2NAL
        outDir = self.home / Path('PAL2NAL')
        Path.mkdir(outDir, exist_ok=True)
        output_file = str(outDir / Path(output_file))

        # Create an alignment for paml input
        P2Ncmd = Pal2NalCommandline(pepaln=aa_alignment, nucfasta=na_fasta, output_file=output_file + '.paml.aln',
                                    output='paml', nogap=nogap, nomismatch=nomismatch)
        print(P2Ncmd)
        pal2nal_flag = True
        while pal2nal_flag is True:
            pal2nal = subprocess.Popen([str(P2Ncmd)], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, encoding='utf-8')
            error = pal2nal.stderr.readlines()
            out = pal2nal.stdout.readlines()
            pal2nal.wait()
            if 'ERROR: inconsistency between the following pep and nuc seqs' in error[0]:
                print('Caught the pal2nal error!')
                print(error[0])
                for err in error:
                    if '>' in err:
                        removed.append(err.strip('>' '\n'))
                multi_fasta_manipulator(na_fasta, removed, na_fasta)
                multi_fasta_manipulator(aa_alignment, removed, aa_alignment)
            else:
                pal2nal_flag = False

            print('Error: ' + str(error))
            print('Out: ' + str(out))
        # Create an alignment for iqtree input
        P2Ncmd = Pal2NalCommandline(pepaln=aa_alignment, nucfasta=na_fasta, output_file=output_file + '.iqtree.aln',
                                    output='fasta', nogap=nogap, nomismatch=nomismatch)
        print(P2Ncmd)
        pal2nal = subprocess.Popen([str(P2Ncmd)], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, encoding='utf-8')
        error = pal2nal.stderr.read()
        out = pal2nal.stdout.read()
        pal2nal.wait()

        print('Error: ' + str(error))
        print('Out: ' + str(out))

        print('Align the nucleic acids using the amino acid alignment.')
        print()




