# -*- coding: utf-8 -*-
"""
Clustal Omega is a multiple sequence alignment program that uses seeded
guide trees and HMM profile-profile techniques to generate alignments between
three or more sequences.

"""
# Import the Clustal Omega wrapper from Biopython
from Bio.Align.Applications import ClustalOmegaCommandline

class ClustalO(object):
    """
    This class aligns genes using parameters similar to the default
    parameters. These parameters include 2 additional iterations for the hmm.
    """

    def __init__(self, infile, outfile, logpath, seqtype="DNA", infmt="fasta",
                  outfmt="fasta", clustalpath='clustalo'):

        # Run clustal omega using the multifasta file
        clustalo_cline = ClustalOmegaCommandline(infile=infile,
                                                 cmd=clustalpath,
                                                 outfile=outfile,
                                                 seqtype=seqtype, # "RNA"
                                                 max_hmm_iterations=2, # Notable
                                                 infmt=infmt,
                                                 outfmt=outfmt, # "aln", "phy"
                                                 iterations=3,  # Notable
                                                 verbose=True,
                                                 force=True,
                                                 log=logpath)
        stdout, stderr = clustalo_cline()
        clustalo_cline()
