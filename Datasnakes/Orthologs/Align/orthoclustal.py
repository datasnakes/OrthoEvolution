"""Clustal Omega is a multiple sequence alignment program that uses seeded
guide trees and HMM profile-profile techniques to generate alignments between
three or more sequences.
"""
# Import the Clustal Omega wrapper from Biopython
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
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
                                                 seqtype="DNA",  # "RNA"
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


class AlignmentEditor:
    """Edit alignments and ensure divisibility by 3."""
    def __init__(self, seqfile, seqtype):
        """Initialize the record handle."""
        self.seqfile = seqfile
        self.seqtype = seqtype
        records = list(SeqIO.parse(self.seqfile, self.seqtype))
        self.records = records

    def divby3(self):
        for record in self.records:
            if len(record.seq) % 3 != 0:
                print('Sequences in %s are not divisible by 3.' % self.seqfile)
                break
            else:
                print('Sequences in %s are divisible by 3.' % self.seqfile)
                break
