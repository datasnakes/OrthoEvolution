"""This is the test suite for Orthologs."""
import unittest
from shutil import rmtree

from OrthoEvol.Orthologs.Blast import OrthoBlastN


class TestOrthologs(unittest.TestCase):
    """Test the Orthologs module."""

    def setUp(self, project="orthology-gpcr"):
        self.project = project

    def delete_project(self, project_path):
        rmtree(project_path)

    def test_orthoblastn(self):
        """Test the OrthoBlastN class."""
        with self.assertRaises(EnvironmentError):
            gpcr_blastn = OrthoBlastN(project=self.project, method=3,
                                      save_data=True, acc_file="gpcr.csv",
                                      copy_from_package=True)
            self.assertEqual(gpcr_blastn.acc_file, "gpcr.csv")
            self.assertTrue(gpcr_blastn.copy_from_package)
        self.delete_project(project_path=self.project)


if __name__ == '__main__':
    unittest.main()
