"""This is the test suite for Orthologs."""
import unittest
from shutil import rmtree
import os

from OrthoEvol.Orthologs.Blast import BaseBlastN, OrthoBlastN
from OrthoEvol.Orthologs.Phylogenetics.PhyML import PhyML
from OrthoEvol.Orthologs.Phylogenetics.TreeViz import TreeViz

class TestOrthologs(unittest.TestCase):
    """Test the Orthologs module."""

    def setUp(self, project="gpcr", project_path="projects"):
        self.project = project
        self.project_path = project_path
        self.cur_dir = os.path.dirname(os.path.abspath(__file__))
        self.join = os.path.join

    def delete_project(self, project_path):
        rmtree(project_path)

    def delete_phyml_output(self):
        os.remove(self.join(self.cur_dir, 'test_data/test.phy_phyml_stats.txt'))
        os.remove(self.join(self.cur_dir, 'test_data/test.phy_phyml_tree.txt'))

    def delete_treeviz_output(self):
        os.remove('example.png')

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
        p = PhyML(infile=self.join(self.cur_dir, 'test_data/test.phy'), datatype='nt')
        p.run(bootstrap=0)
        self.assertIsNotNone(self.join(self.cur_dir, 'test_data/test.phy_phyml_stats.txt'))
        self.assertIsNotNone(self.join(self.cur_dir, 'test_data/test.phy_phyml_tree.txt'))
        self.delete_phyml_output()

    def test_treeviz(self):
        """Test the TreeViz class."""
        t = TreeViz(path=self.join(self.cur_dir, 'test_data/test_tree.txt'),
                    tree_format='newick')
        t.draw_tree()
        t.save_tree('example.png')
        self.assertIsNotNone('example.png')
        self.delete_treeviz_output()
        
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
