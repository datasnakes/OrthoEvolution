"""Test the PhyML executable.
https://github.com/biopython/biopython/blob/master/Tests/test_phyml_tool.py
"""
import sys
import os
import unittest
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio import MissingExternalDependencyError


class PhymlTest(unittest.TestCase):
    """Tests for application wrappers."""

    def __init__(self):
        # Try to avoid problems when the OS is in another language
        os.environ['LANG'] = 'C'

        phyml_exe = None
        exe_name = "PhyML-3.1_win32.exe" if sys.platform == "win32" else "phyml"
        from Bio._py3k import getoutput
        try:
            output = getoutput(exe_name + " --version")
            if "not found" not in output and "20" in output:
                phyml_exe = exe_name
        except OSError:
            # Python 2.6 or 2.7 on Windows XP:
            # WindowsError: [Error 2] The system cannot find the file specified
            # Python 3.3 or 3.4 on Windows XP:
            # FileNotFoundError: [WinError 2] The system cannot find the file
            # specified
            pass

        if not phyml_exe:
            raise MissingExternalDependencyError(
                "Install PhyML 3.0 if you want to use the \
                Bio.Phylo.Applications wrapper.")

        # Example Phylip file with 13 aligned protein sequences
        EX_PHYLIP = 'HTR1E_aligned.phy'
        return EX_PHYLIP

    def test_phyml(self):
        """Run PhyML using the wrapper."""
        cmd = PhymlCommandline(
            self.phyml_exe,
            input=self.EX_PHYLIP,
            datatype='nt')
        # Smoke test
        try:
            out, err = cmd()
            self.assertTrue(len(out) > 0)
            self.assertEqual(len(err), 0)
            # Check the output tree
            tree = Phylo.read(self.EX_PHYLIP + '_phyml_tree.txt', 'newick')
            self.assertEqual(tree.count_terminals(), 13)
        finally:
            # Clean up generated files
            for suffix in ['_phyml_tree.txt', '_phyml_stats.txt']:
                fname = self.EX_PHYLIP + suffix
                if os.path.isfile(fname):
                    os.remove(fname)
