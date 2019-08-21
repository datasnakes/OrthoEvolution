import sys
import shutil

from Bio.Phylo.Applications import PhymlCommandline
from Bio.Application import ApplicationError
from Bio import AlignIO

from OrthoEvol.Tools.logit import LogIt


class PhyML(object):
    """The PhyML class uses Biopython's PhyMLCommandline wrapper to generate
    trees from the PhyML executable."""

    def __init__(self, infile, datatype="aa"):
        """Input a phylip formatted alignment file and specify a datatype.

        :param infile: An input file that is phylip formatted.
        :type infile: str
        :param datatype: The datatype of the infile ("nt"/"aa"), defaults to "aa"
        :type datatype: str, optional

        If you're using Linux, ensure that your phyml path is set in your bash
        profile. If you're using Windows, this function will look for the name
        of the executable 'PhyML-3.1_win32.exe'.
        """
        # Set up logging
        self.phyml_log = LogIt().default(logname="Phyml", logfile=None)
        # Check that the phyml executable is in the path
        self.phyml_exe = self._check_exe()
        self.datatype = datatype
        if self._validate_format(infile):
            self.infile = infile

    def _validate_format(self):
        """"Validate the format of the input file.

        :param infile: An input file that is phylip formatted.
        :type infile: str
        """
        try:
            AlignIO.read(open(self.infile), "phylip")
        except ValueError as e:
            self.phyml_log.exception(e)
        else:
            return True
        return False

    def _check_exe(self):
        """Check to see if the phyml exe is in the path."""
        phyml_exe = None
        win32 = "win32"
        executable = "PhyML-3.1_win32.exe"
        exe_name = executable if sys.platform == win32 else "phyml"
        phyml_exe = exe_name
        if shutil.which(phyml_exe):
            return phyml_exe
        else:
            self.phyml_log.error("%s is not in the path." % phyml_exe)

    def run(self, model="WAG", alpha="e", bootstrap=100):
        """"Run phyml."""
        try:
            run_phyml = PhymlCommandline(self.phyml_exe,
                                         input=self.infile,
                                         datatype=self.datatype, model=model,
                                         alpha=alpha, bootstrap=bootstrap)
            self.phyml_log.info("Running %s on %s" % (self.phyml_exe,
                                                      self.infile))
            out_log, err_log = run_phyml()
            if out_log:
                self.phyml_log.info(out_log)
            if err_log:
                self.phyml_log.error(err_log)
        except ApplicationError as e:
            self.phyml_log.exception(e)
