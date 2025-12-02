import sys
import shutil

from Bio.Phylo.Applications import PhymlCommandline
try:
    from Bio.Application import ApplicationError
except ImportError:
    # Bio.Application is deprecated in newer biopython versions
    # Use subprocess.CalledProcessError as fallback
    from subprocess import CalledProcessError as ApplicationError
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
        # Validate format and set infile
        if not self._validate_format(infile):
            raise ValueError(f"Invalid phylip format for file: {infile}")
        self.infile = infile

    def _validate_format(self, infile):
        """Validate the format of the input file.

        :param infile: An input file that is phylip formatted.
        :type infile: str
        """
        try:
            with open(infile, 'r') as f:
                AlignIO.read(f, "phylip")
            return True
        except (ValueError, IOError) as e:
            self.phyml_log.exception(e)
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
            error_message = f"{phyml_exe} is not in the PATH. Please ensure that the PhyML executable is installed and available in your system's PATH."
            self.phyml_log.error(error_message)
            raise FileNotFoundError(error_message)

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
