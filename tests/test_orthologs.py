"""This is the test suite for Orthologs."""
import unittest
from shutil import rmtree

from OrthoEvol.Orthologs.Blast import BaseBlastN


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
            gpcr_blastn = BaseBlastN(project=self.project, method=3,
                                     save_data=True, acc_file="gpcr.csv",
                                     proj_mana=None,
                                     project_path=self.project_path,
                                     copy_from_package=True)
            self.assertEqual(gpcr_blastn.proj_mana, None)
            self.assertEqual(gpcr_blastn.acc_file, "gpcr.csv")
            self.assertTrue(gpcr_blastn.copy_from_package)
            self.delete_project(project_path=self.project_path)


if __name__ == '__main__':
    unittest.main()
