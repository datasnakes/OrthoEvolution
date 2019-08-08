"""This is the test suite for Orthologs."""
import unittest
from shutil import rmtree
import os

from OrthoEvol.Orthologs.Blast import BaseBlastN
from OrthoEvol.Orthologs.Phylogenetics.PhyML import PhyML


class TestOrthologs(unittest.TestCase):
    """Test the Orthologs module."""

    def setUp(self, project="gpcr", project_path="projects"):
        self.project = project
        self.project_path = project_path

    def delete_project(self, project_path):
        rmtree(project_path)

    def delete_phyml_output(self):
        os.remove('test_data/HTR1E_aligned.phy_phyml_stats.txt')
        os.remove('test_data/HTR1E_aligned.phy_phyml_tree.txt')

    def test_baseblastn(self):
        """Test the BaseBlastN class."""
        # The with statement is for travisci where a BLASTDB variable
        # is not set.
        # TIP: Remove the with statement if testing with BLASTDB in your
        # environment variables.
        with self.assertRaises(EnvironmentError):
            gpcr_blastn = BaseBlastN(project=self.project, method=1,
                                     save_data=True, acc_file="gpcr.csv",
                                     copy_from_package=True,
                                     ref_species='Homo_sapiens',
                                     proj_mana=None,
                                     project_path=self.project_path)
            self.assertEqual(gpcr_blastn.proj_mana, None)
            self.assertEqual(gpcr_blastn.acc_file, "gpcr.csv")
            self.assertTrue(gpcr_blastn.copy_from_package)
            self.delete_project(project_path=self.project_path)

    def test_phyml(self):
        """Test the PhyML class."""
        PhyML(infile='test_data/HTR1E_aligned.phy', datatype='aa').run()
        self.assertIsNotNone('test_data/HTR1E_aligned.phy_phyml_stats.txt')
        self.assertIsNotNone('test_data/HTR1E_aligned.phy_phyml_tree.txt')
        self.delete_phyml_output()


if __name__ == '__main__':
    unittest.main()
