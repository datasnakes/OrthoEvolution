import sys

from Bio.Phylo.Applications import PhymlCommandline

from OrthoEvol.Tools.logit import LogIt


class PhyML(object):
    """The PhyML class uses Biopython's PhyMLCommandline wrapper to generate 
    trees from the PhyML executable."""

    def __init__(self, phyml_input, datatype="aa"):
        """Run phyml to generate tree results.

        If you're using Linux, ensure that your phyml path is set in your bash
        profile. If you're using Windows, this function will look for the name
        of the executable 'PhyML-3.1_win32.exe'.
        """
        self.phyml_log = LogIt().default(logname="GenBank", logfile=None)

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
        self.phyml_input = phyml_input

    def run(self):
        """"Run phyml.

        Input a phylip formatted alignment file and describe the datatype
        ('nt' or 'aa').
        """

        run_phyml = PhymlCommandline(self.phyml_exe,
                                     input=self.phyml_input,
                                     datatype=self.datatype)
        out_log, err_log = run_phyml()
        self.phyml_log(out_log)
        self.phyml_log(err_log)
