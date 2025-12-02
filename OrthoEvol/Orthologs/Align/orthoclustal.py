"""Clustal Omega is a multiple sequence alignment program that uses seeded
guide trees and HMM profile-profile techniques to generate alignments between
three or more sequences.
"""
from Bio.Align.Applications import ClustalOmegaCommandline
try:
    from Bio.Application import ApplicationError
except ImportError:
    # Bio.Application is deprecated in newer biopython versions
    # Use subprocess.CalledProcessError as fallback
    from subprocess import CalledProcessError as ApplicationError

from OrthoEvol.Tools.logit import LogIt

stop_codons = ['TAG', 'TAA', 'TGA']


class ClustalO(object):
    """Align genes using Clustal Omega.

    This class is a further wrapper around Biopython's ClustalOmegaCommandline.

    :param infile: Path/Name of multiple fasta file.
    :type infile: str
    :param outfile: Path/Name of multiple alignment file.
    :type outfile: str
    :param logpath: Path to logfile.
    :type logpath: str or None
    :param outfmt: Format of the output multiple alignment file (e.g., 'fasta', 'clustal', 'phylip').
    :type outfmt: str
    """

    def __init__(self, infile, outfile, logpath=None, outfmt="fasta"):
        """Set up the logger and the parameters.

        :param infile: Path/Name of multiple fasta file.
        :type infile: str
        :param outfile: Path/Name of multiple alignment file.
        :type outfile: str
        :param logpath: Path to logfile.
        :type logpath: str or None
        :param outfmt: Format of the output multiple alignment file.
        :type outfmt: str
        """
        self.infile = infile
        self.outfile = outfile
        self.outfmt = outfmt
        self.logpath = logpath
        self.clustalolog = LogIt().default('clustalo', logfile=self.logpath)

    def runclustalomega(self):
        """Run Clustal Omega alignment.

        Executes the Clustal Omega command line tool to perform multiple
        sequence alignment on the input file.
        """

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
