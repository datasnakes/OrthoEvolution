"""Clustal Omega is a multiple sequence alignment program that uses seeded
guide trees and HMM profile-profile techniques to generate alignments between
three or more sequences.
"""
# Import the Clustal Omega wrapper from Biopython
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio import AlignIO
stop_codons = ['TAG', 'TAA', 'TGA']


class ClustalO:
    """This class aligns genes using parameters similar to the default parameters.

    These parameters include 2 additional iterations for the hmm.
    """

    def __init__(self, infile, outfile, logpath, outfmt="fasta"):

        # Run clustal omega using the multifasta file
        clustalo_cline = ClustalOmegaCommandline(infile=infile,
                                                 cmd="clustalo",
                                                 outfile=outfile,
                                                 seqtype="PROTEIN",  # "RNA"/"DNA"
                                                 max_hmm_iterations=2,
                                                 infmt="fasta",
                                                 outfmt=outfmt,  # "aln", "phy"
                                                 iterations=3,  # Notable
                                                 verbose=True,
                                                 force=True,
                                                 log=logpath)
        stdout, stderr = clustalo_cline()
        clustalo_cline()
        if stderr:
            print(stderr)
        if stdout:
            print(stdout)


#class SequenceEditor:
#    """Edit sequences prior to running Clustal Omega."""
#    def __init__(self, seqfile, seqtype):
#        """Initialize the sequence file record."""
#        self.seqfile = seqfile
#        self.seqtype = seqtype
#        records = list(SeqIO.parse(self.seqfile, self.seqtype))
#        self.records = records
#
#
#class AlignmentEditor:
#    """Edit alignments to ensure length is a multiple of 3 for PAML."""
#    def __init__(self, seqfile, seqtype):
#        """Initialize the multiple sequence alignment file record."""
#        self.seqfile = seqfile
#        self.seqtype = seqtype
#        self.alignments = AlignIO.parse(self.seqfile, self.seqtype)
#
#    def divby3(self):
#        """Ensure divisibility by 3."""
#        for alignment in self.alignments:
#            if alignment.get_alignment_length() % 3 != 0:
#                print('Sequences in %s are not divisible by 3.' % alignment)
#                break
#
#    def pamlslice(self):
#        """Slice or add to alignments to ensure multiple of 3 for PAML."""
#        print(self.pamlslice.__doc__)
