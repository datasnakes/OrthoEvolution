from Datasnakes.Orthologs.Align.QualityControl.guidance2 import Guidance2Commandline
from Datasnakes.Orthologs.Align.QualityControl.pal2nal import Pal2NalCommandline
from Datasnakes.Orthologs.Genbank.genbank import GenBank
from Datasnakes.Orthologs.Genbank.utils import multi_fasta_manipulator
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from pathlib import Path
import os
import subprocess


class Alignment(GenBank):

    def __init__(self, program, **kwargs):
        super().__init__(**kwargs)


        print()

    def guidance2(self, seqFile, msaProgram, seqType, outDir="GUIDANCE2", seqFilter=None, columnFilter=None, maskFilter=None, **kwargs):
        # Name and Create the output directory
        gene = Path(seqFile).stem
        geneDir = self.raw_data / Path(gene)
        outDir = self.raw_data / Path(gene) / Path(outDir)
        Path.mkdir(outDir, parents=True, exist_ok=True)

        # Add the Guidance 2 cutoffs to the keyword arguments
        if 'seqCutoff' not in kwargs.keys():
            kwargs['seqCutoff'] = 0.6
        if 'colCutoff' not in kwargs.keys():
            kwargs['colCutoff'] = 0.93

        # filter = "inclusive", "exclusive", "neutral", or None (default) or "mask" (columnFilter)
        if seqFilter is not None:

            # Filter "bad" sequences and iterate over the good ones until the cutoff is reached.
            iter_flag = True
            iteration = 0
            while iter_flag is True:

                iteration += 1
                # Create paths for output files
                iterDir = Path(outDir) / Path('iter_%s' % iteration)
                g2_rem_file = str(iterDir / Path('Seqs.Orig.fas.FIXED.Removed_Seq.With_Names'))  # Need for all iterations
                g2_seqFile = str(geneDir / Path(self.gene + '_G2.ffn'))  # Need for all iterations
                rem_file = str(geneDir / Path(self.gene + '_G2_removed.ffn'))   # Need for all iterations
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
                        # Filter the input NA fasta file using the updated rem_file
                        multi_fasta_manipulator(seqFile, rem_file, g2_seqFile, manipulation='remove')
                        iterFlag = True
                    else:
                        iterFlag = False

        # Filter the "bad" columns using good sequences
        if columnFilter is not None:
            print()
        # Mask the "bad" residues using the good sequences
        elif maskFilter is not None:
            print()





