from Bio.Phylo.Applications import PhymlCommandline
import sys


class PhyML(object):
    """
    The PhyML class uses Biopython's PhyMLCommandline wrapper to generate trees
    from the PhyML executable. This class also converts a fasta formatted
    multiple sequence alignment file into relaxed phylip format.
    """

    def __init__(self, phyml_input, datatype='nt'):
        """Run phyml to generate tree results.

        If you're using Linux, ensure that your phyml path is set in your bash
        profile. If you're using Windows, this function will look for the name
        of the executable 'PhyML-3.1_win32.exe'.
        """
        # Use the phyml executable file
        phyml_exe = None

        # This is mainly intended for windows use or use with an executable
        # file
        exe_name = "PhyML-3.1_win32.exe" if sys.platform == "win32" else "phyml"
        phyml_exe = exe_name
        phyml_exe = phyml_exe

        # Create the command & run phyml
        # Input a phylip formatted alignment file and describe the datatype
        # ('nt' or 'aa')
        run_phyml = PhymlCommandline(
            phyml_exe, input=phyml_input, datatype=datatype)
        print(run_phyml)
        out_log, err_log = run_phyml()
