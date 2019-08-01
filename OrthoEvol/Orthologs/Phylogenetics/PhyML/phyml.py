import sys

from Bio.Phylo.Applications import PhymlCommandline

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
        self.phyml_log = LogIt().default(logname="Phyml", logfile=None)

        # Use the phyml executable file
        phyml_exe = None

        # This is mainly intended for windows use or use with an executable
        # file
        win32 = "win32"
        executable = "PhyML-3.1_win32.exe"
        exe_name = executable if sys.platform == win32 else "phyml"
        phyml_exe = exe_name
        self.phyml_exe = phyml_exe
        self.datatype = datatype
        self.phyml_input = infile

    def _validate_format(self, infile):
        """"Validate the format of the input file.

        :param infile: An input file that is phylip formatted.
        :type infile: str
        """
        pass

    def _check_exe(self):
        """Check to see if the phyml exe is in the path."""
        pass

    def run(self):
        """"Run phyml."""
        # TODO: Add try/except logic.
        run_phyml = PhymlCommandline(self.phyml_exe,
                                     input=self.phyml_input,
                                     datatype=self.datatype)
        out_log, err_log = run_phyml()
        if out_log:
            self.phyml_log.info(out_log)
        if err_log:
            self.phyml_log.error(err_log)
