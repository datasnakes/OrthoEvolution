"""Clustal Omega is a multiple sequence alignment program that uses seeded
guide trees and HMM profile-profile techniques to generate alignments between
three or more sequences.
"""
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Application import ApplicationError

from Datasnakes.Tools import LogIt

stop_codons = ['TAG', 'TAA', 'TGA']


class ClustalO(object):
    """Align genes using Clustal Omega.

    This class is a further wrapper around Biopython's ClustalOmegaCommandline.

    :param infile:
    :param outfile:
    :param logpath:
    :param outfmt:
    """

    def __init__(self, infile, outfile, logpath=None, outfmt="fasta"):
        """Set up the logger and the parameters."""
        self.infile = infile
        self.outfile = outfile
        self.outfmt = outfmt
        self.logpath = logpath
        self.clustalolog = LogIt().default('clustalo', logfile=self.logpath)

    def runclustalomega(self):
        """Run clustalomega."""
        try:
            # Run clustal omega using the multifasta file
            clustalo_cline = ClustalOmegaCommandline(infile=self.infile,
                                                     cmd="clustalo",
                                                     outfile=self.outfile,
                                                     # "RNA"/"DNA"
                                                     seqtype="PROTEIN",
                                                     max_hmm_iterations=2,
                                                     infmt="fasta",
                                                     # "aln", "phy"
                                                     outfmt=self.outfmt,
                                                     iterations=3,  # Notable
                                                     verbose=True,
                                                     force=True,
                                                     log=self.logpath)
            clustalo_cline()
            stdout, _ = clustalo_cline()
            self.clustalolog.info(stdout)

        except ApplicationError as err:
            self.clustalolog.error(err)
