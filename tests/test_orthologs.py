"""This is the test suite for Orthologs."""
import unittest
from shutil import rmtree

from OrthoEvol.Orthologs.Blast import BaseBlastN, OrthoBlastN


class TestOrthologs(unittest.TestCase):
    """Test the Orthologs module."""

    def setUp(self, project="gpcr", project_path="projects"):
        self.project = project
        self.project_path = project_path

    def delete_project(self, project_path):
        rmtree(project_path)

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

    def test_orthoblastn(self):
        """Test the OrthoBlastN class."""
        with self.assertRaises(EnvironmentError):
            ortho_blastn = OrthoBlastN(project="orthology-project",
                                       method=1, save_data=True,
                                       acc_file="gpcr.csv",
                                       copy_from_package=True)
            self.assertEqual(ortho_blastn.ref_species, 'Homo_sapiens')
            self.assertTrue(ortho_blastn.copy_from_package)


if __name__ == '__main__':
    unittest.main()
