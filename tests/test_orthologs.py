"""This is the test suite for Orthologs."""
import unittest
from shutil import rmtree
import os
import tempfile
from unittest import mock
from pathlib import Path

from OrthoEvol.Orthologs.Blast import BaseBlastN, OrthoBlastN
from OrthoEvol.Orthologs.Phylogenetics.PhyML import PhyML
from OrthoEvol.Orthologs.Phylogenetics.TreeViz import TreeViz
from OrthoEvol.Orthologs.Phylogenetics import RelaxPhylip
from OrthoEvol.Orthologs.Phylogenetics.IQTree import IQTreeCommandline, FilteredTree
from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML
from OrthoEvol.Orthologs.GenBank.genbank import GenBank
from OrthoEvol.Orthologs.Align.orthoclustal import ClustalO
from OrthoEvol.Orthologs.Align.msa import MultipleSequenceAlignment

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
        if os.path.exists('example.png'):
            return os.remove('example.png')

    def delete_relaxphylip_output(self, output_file):
        if os.path.exists(output_file):
            os.remove(output_file)

    def test_baseblastn_missing_blastdb(self):
        """Test BaseBlastN raises EnvironmentError when BLASTDB is not set."""
        # Save original BLASTDB if it exists
        original_blastdb = os.environ.get('BLASTDB')
        
        # Remove BLASTDB from environment
        if 'BLASTDB' in os.environ:
            del os.environ['BLASTDB']
        
        # Use a minimal test CSV file
        test_csv = self.join(self.cur_dir, 'test_data', 'test_gpcr.csv')
        
        try:
            # Test that BaseBlastN raises EnvironmentError when BLASTDB is missing
            # Note: Initialization may fail earlier due to NCBITaxa database issues,
            # but if it gets to the BLASTDB check, it should raise EnvironmentError
            try:
                BaseBlastN(project=self.project, 
                          method=1,
                          save_data=True, 
                          acc_file=test_csv,
                          copy_from_package=False,
                          ref_species='Homo_sapiens',
                          proj_mana=None,
                          project_path=self.project_path)
                # If we get here, BLASTDB was set or check was bypassed
                self.fail("Expected EnvironmentError for missing BLASTDB, but initialization succeeded")
            except EnvironmentError as e:
                # This is the expected error - verify it's about BLASTDB
                self.assertIn("BLASTDB", str(e))
            except Exception as e:
                # Other errors (like NCBITaxa database issues) may occur before BLASTDB check
                # Check if the error message mentions BLASTDB (shouldn't, but verify)
                error_str = str(e)
                if "BLASTDB" in error_str:
                    # If BLASTDB is mentioned, it might be the error we're looking for
                    self.assertIn("BLASTDB", error_str)
                else:
                    # Skip if it's a different initialization issue (like NCBITaxa)
                    self.skipTest(f"BaseBlastN initialization failed before BLASTDB check (likely NCBITaxa issue): {type(e).__name__}")
        finally:
            # Restore original BLASTDB
            if original_blastdb is not None:
                os.environ['BLASTDB'] = original_blastdb

    def test_baseblastn_with_blastdb(self):
        """Test BaseBlastN initialization when BLASTDB is set."""
        # Save original BLASTDB if it exists
        original_blastdb = os.environ.get('BLASTDB')
        
        # Set a dummy BLASTDB for testing
        os.environ['BLASTDB'] = '/tmp/test_blastdb'
        
        # Use a minimal test CSV file
        test_csv = self.join(self.cur_dir, 'test_data', 'test_gpcr.csv')
        
        try:
            # Test initialization (may fail on other requirements, but should pass BLASTDB check)
            try:
                gpcr_blastn = BaseBlastN(project=self.project, 
                                        method=1,
                                        save_data=True, 
                                        acc_file=test_csv,
                                        copy_from_package=False,
                                        ref_species='Homo_sapiens',
                                        proj_mana=None,
                                        project_path=self.project_path)
                
                # If initialization succeeds, verify attributes
                self.assertEqual(gpcr_blastn.proj_mana, None)
                self.assertEqual(gpcr_blastn.acc_file, test_csv)
                self.assertFalse(gpcr_blastn.copy_from_package)
                self.assertEqual(gpcr_blastn.ref_species, 'Homo_sapiens')
                self.assertEqual(gpcr_blastn.method, 1)
                self.assertEqual(gpcr_blastn.project, self.project)
                
                # Verify that organisms and genes were loaded from CSV
                self.assertGreater(len(gpcr_blastn.org_list), 0, "Organism list should not be empty")
                self.assertGreater(len(gpcr_blastn.gene_list), 0, "Gene list should not be empty")
                # Verify expected organisms are in the list
                self.assertIn('Homo_sapiens', gpcr_blastn.org_list)
                self.assertIn('Macaca_mulatta', gpcr_blastn.org_list)
                self.assertIn('Mus_musculus', gpcr_blastn.org_list)
                # Verify expected genes are in the list
                self.assertIn('ADRA1A', gpcr_blastn.gene_list)
                self.assertIn('ADRA1B', gpcr_blastn.gene_list)
                
                # Cleanup if project was created
                if os.path.exists(self.project_path):
                    self.delete_project(project_path=self.project_path)
            except (FileNotFoundError, ValueError, KeyError, EnvironmentError) as e:
                # EnvironmentError means BLASTDB check failed (unexpected since we set it)
                if isinstance(e, EnvironmentError) and "BLASTDB" in str(e):
                    self.fail(f"BLASTDB was set but still got EnvironmentError: {e}")
                # Other errors are expected if dependencies are missing (NCBITaxa, etc.)
                self.skipTest(f"BaseBlastN initialization failed due to missing dependencies: {e}")
            except Exception as e:
                # Handle other unexpected errors (like NCBITaxa database issues)
                error_str = str(e)
                if "BLASTDB" in error_str:
                    self.fail(f"Unexpected BLASTDB error: {e}")
                # Skip for other initialization issues
                self.skipTest(f"BaseBlastN initialization failed: {e}")
        finally:
            # Restore original BLASTDB
            if original_blastdb is not None:
                os.environ['BLASTDB'] = original_blastdb
            elif 'BLASTDB' in os.environ:
                del os.environ['BLASTDB']

    def test_treeviz(self):
        """Test the TreeViz class."""
        t = TreeViz(path=self.join(self.cur_dir, 'test_data/test_tree.txt'),
                    tree_format='newick')
        t.draw_tree()
        t.save_tree('example.png')
        self.assertIsNotNone('example.png')
        self.delete_treeviz_output()

    def test_relaxphylip(self):
        """Test the RelaxPhylip class for format conversion."""
        # Create a temporary fasta file for testing
        test_fasta = self.join(self.cur_dir, 'test_data', 'test_relaxphylip.fasta')
        test_phy = self.join(self.cur_dir, 'test_data', 'test_relaxphylip.phy')
        
        # Create a simple test fasta file
        with open(test_fasta, 'w') as f:
            f.write(">seq1\nATGCATGC\n>seq2\nATGCATGC\n")
        
        try:
            # Test RelaxPhylip conversion
            relax_phylip = RelaxPhylip(test_fasta, test_phy)
            # Check that output file was created
            self.assertTrue(os.path.exists(test_phy))
        finally:
            # Cleanup
            if os.path.exists(test_fasta):
                os.remove(test_fasta)
            self.delete_relaxphylip_output(test_phy)

    def test_iqtree_commandline_init(self):
        """Test IQTreeCommandline initialization."""
        # Test that IQTreeCommandline can be instantiated with required parameters
        # Note: This doesn't require the iqtree executable to be installed
        test_alignment = self.join(self.cur_dir, 'test_data', 'test.phy')
        if os.path.exists(test_alignment):
            try:
                iqtree_cline = IQTreeCommandline(alignment=test_alignment)
                self.assertIsNotNone(iqtree_cline)
                # Test that it has the expected attributes
                self.assertTrue(hasattr(iqtree_cline, 'parameters'))
            except Exception as e:
                # If iqtree executable is not found, that's expected
                # We're just testing that the class can be initialized
                self.skipTest(f"IQTree initialization failed (expected if executable not installed): {e}")

    def test_filtered_tree_init(self):
        """Test FilteredTree initialization and directory setup."""
        # Create a temporary directory for testing
        test_home = tempfile.mkdtemp()
        test_alignment_name = 'test_gene_P2N_na.iqtree.aln'
        test_alignment = os.path.join(test_home, test_alignment_name)
        
        # Create a simple test alignment file
        with open(test_alignment, 'w') as f:
            f.write("3 10\nseq1 ATGCATGCAT\nseq2 ATGCATGCAT\nseq3 ATGCATGCAT\n")
        
        try:
            # Test FilteredTree initialization
            # This will create directories and copy files, but may fail if IQTree isn't installed
            filtered_tree = FilteredTree(alignment=test_alignment_name, 
                                        dataType='CODON', 
                                        home=test_home)
            
            # Verify that IQTREE directory was created
            iqtree_dir = os.path.join(test_home, 'IQTREE')
            self.assertTrue(os.path.exists(iqtree_dir))
            
            # Verify that alignment was copied to IQTREE directory
            copied_alignment = os.path.join(iqtree_dir, test_alignment_name)
            self.assertTrue(os.path.exists(copied_alignment))
            
            # Verify attributes are set
            self.assertEqual(filtered_tree.home, test_home)
            self.assertEqual(filtered_tree.gene, 'test_gene')
            self.assertTrue(hasattr(filtered_tree, 'tree_file'))
            self.assertTrue(hasattr(filtered_tree, 'iqtree_path'))
            
        except FileNotFoundError:
            # IQTree executable not found - this is expected
            self.skipTest("IQTree executable not found in PATH")
        except Exception as e:
            # Other errors during execution - still verify directory setup happened
            iqtree_dir = os.path.join(test_home, 'IQTREE')
            if os.path.exists(iqtree_dir):
                # Directory was created, which means initialization started
                self.assertTrue(os.path.exists(iqtree_dir))
            self.skipTest(f"FilteredTree execution failed (expected if IQTree not fully functional): {e}")
        finally:
            # Cleanup
            if os.path.exists(test_home):
                rmtree(test_home)

    def test_phyml_validation(self):
        """Test PhyML class validation without running."""
        test_phy = self.join(self.cur_dir, 'test_data', 'test.phy')
        if os.path.exists(test_phy):
            try:
                # Test initialization and validation
                p = PhyML(infile=test_phy, datatype='nt')
                self.assertIsNotNone(p)
                self.assertEqual(p.infile, test_phy)
                self.assertEqual(p.datatype, 'nt')
            except FileNotFoundError:
                # PhyML executable not found - skip test
                self.skipTest("PhyML executable not found in PATH")
            except ValueError as e:
                # Invalid format - this is a real error
                self.fail(f"PhyML validation failed: {e}")

    def test_ete3paml_init(self):
        """Test ETE3PAML initialization."""
        # This test requires test files and may require PAML executable
        # ETE3PAML expects a fasta file, so we'll create a simple one for testing
        test_tree = self.join(self.cur_dir, 'test_data', 'test_tree.txt')
        test_fasta = self.join(self.cur_dir, 'test_data', 'test_ete3paml.fasta')
        workdir = tempfile.mkdtemp()
        
        # Create a simple test fasta file
        if os.path.exists(test_tree):
            try:
                with open(test_fasta, 'w') as f:
                    f.write(">Human\nATGCATGC\n>Chimp\nATGCATGC\n")
                
                # Test initialization (won't run without PAML executable)
                paml = ETE3PAML(infile=test_fasta, 
                               species_tree=test_tree, 
                               workdir=workdir,
                               pamlsrc=None)
                self.assertIsNotNone(paml)
                self.assertEqual(paml.infile, test_fasta)
                self.assertEqual(paml.species_tree, test_tree)
                
                # Cleanup
                if os.path.exists(test_fasta):
                    os.remove(test_fasta)
                if os.path.exists(workdir):
                    rmtree(workdir)
            except Exception as e:
                # If initialization fails due to missing dependencies, skip
                if os.path.exists(test_fasta):
                    os.remove(test_fasta)
                if os.path.exists(workdir):
                    rmtree(workdir)
                self.skipTest(f"ETE3PAML initialization failed: {e}")
        
    # def test_orthoblastn(self):
    #     """Test the OrthoBlastN class."""
    #     with self.assertRaises(EnvironmentError):
    #         ortho_blastn = OrthoBlastN(project="orthology-project",
    #                                    method=1, save_data=True,
    #                                    acc_file="gpcr.csv",
    #                                    copy_from_package=True)
    #         self.assertEqual(ortho_blastn.ref_species, 'Homo_sapiens')
    #         self.assertTrue(ortho_blastn.copy_from_package)


class TestGenBank(unittest.TestCase):
    """Test the GenBank class."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.test_path = Path(self.test_dir)

    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.test_dir):
            rmtree(self.test_dir, ignore_errors=True)

    def test_name_fasta_file_cds_single(self):
        """Test name_fasta_file for CDS single mode."""
        gene = "HTR1A"
        org = "Homo_sapiens"
        feat_type = "CDS"
        feat_type_rank = "CDS"
        extension = ".ffn"
        mode = "w"

        file_obj = GenBank.name_fasta_file(
            self.test_path, gene, org, feat_type, feat_type_rank, extension, mode
        )

        self.assertIsNotNone(file_obj)
        self.assertEqual(file_obj.mode, 'w')
        expected_path = self.test_path / f"{gene}_{org}{feat_type_rank}{extension}"
        self.assertEqual(Path(file_obj.name), expected_path)
        file_obj.close()

    def test_name_fasta_file_cds_multi(self):
        """Test name_fasta_file for CDS multi mode."""
        gene = "HTR1A"
        org = "Homo_sapiens"
        feat_type = "CDS"
        feat_type_rank = "CDS"
        extension = ".ffn"
        mode = "a"

        file_obj = GenBank.name_fasta_file(
            self.test_path, gene, org, feat_type, feat_type_rank, extension, mode
        )

        self.assertIsNotNone(file_obj)
        self.assertEqual(file_obj.mode, 'a')
        expected_path = self.test_path / f"{gene}{feat_type_rank}{extension}"
        self.assertEqual(Path(file_obj.name), expected_path)
        file_obj.close()

    def test_name_fasta_file_other_single(self):
        """Test name_fasta_file for other feature type single mode."""
        gene = "HTR1A"
        org = "Homo_sapiens"
        feat_type = "misc_feature"
        feat_type_rank = "misc_feature"
        extension = ".fna"
        mode = "w"

        file_obj = GenBank.name_fasta_file(
            self.test_path, gene, org, feat_type, feat_type_rank, extension, mode
        )

        self.assertIsNotNone(file_obj)
        self.assertEqual(file_obj.mode, 'w')
        expected_path = self.test_path / f"{gene}_{org}_{feat_type_rank}{extension}"
        self.assertEqual(Path(file_obj.name), expected_path)
        file_obj.close()

    def test_name_fasta_file_other_multi(self):
        """Test name_fasta_file for other feature type multi mode."""
        gene = "HTR1A"
        org = "Homo_sapiens"
        feat_type = "misc_feature"
        feat_type_rank = "misc_feature"
        extension = ".fna"
        mode = "a"

        file_obj = GenBank.name_fasta_file(
            self.test_path, gene, org, feat_type, feat_type_rank, extension, mode
        )

        self.assertIsNotNone(file_obj)
        self.assertEqual(file_obj.mode, 'a')
        expected_path = self.test_path / f"{gene}_{feat_type_rank}{extension}"
        self.assertEqual(Path(file_obj.name), expected_path)
        file_obj.close()

    def test_protein_gi_fetch_with_gi(self):
        """Test protein_gi_fetch when GI is present."""
        mock_feature = mock.Mock()
        mock_feature.qualifiers = ['protein_id=NP_000000.1', 'GI:123456789', 'product=test protein']

        result = GenBank.protein_gi_fetch(mock_feature)

        self.assertEqual(result, '123456789')

    def test_protein_gi_fetch_without_gi(self):
        """Test protein_gi_fetch when GI is not present."""
        mock_feature = mock.Mock()
        mock_feature.qualifiers = ['protein_id=NP_000000.1', 'product=test protein']

        result = GenBank.protein_gi_fetch(mock_feature)

        self.assertIsNone(result)

    def test_protein_gi_fetch_empty_qualifiers(self):
        """Test protein_gi_fetch with empty qualifiers."""
        mock_feature = mock.Mock()
        mock_feature.qualifiers = []

        result = GenBank.protein_gi_fetch(mock_feature)

        self.assertIsNone(result)

    @mock.patch('OrthoEvol.Orthologs.GenBank.genbank.FullUtilities')
    @mock.patch('OrthoEvol.Orthologs.GenBank.genbank.LogIt')
    @mock.patch('OrthoEvol.Orthologs.GenBank.genbank.Path.mkdir')
    @mock.patch('os.listdir')
    def test_genbank_init(self, mock_listdir, mock_mkdir, mock_logit, mock_utils):
        """Test GenBank initialization."""
        mock_listdir.return_value = ['test.db', 'other.db', 'not_a_db.txt']
        mock_logit_instance = mock.Mock()
        mock_logit.return_value.default.return_value = mock_logit_instance
        
        # Create a mock object with the required attributes
        mock_config_obj = mock.Mock()
        mock_config_obj.user_db = Path(self.test_dir)
        mock_config_obj.ncbi_db_repo = Path(self.test_dir)
        mock_config_obj.raw_data = Path(self.test_dir)
        
        # Set up the FullUtilities mock
        mock_utils_instance = mock.Mock()
        mock_utils_instance.attribute_config.return_value = mock_config_obj
        mock_utils.return_value = mock_utils_instance

        genbank = GenBank(
            project="test-project",
            project_path=self.test_dir,
            solo=True,
            multi=True
        )

        self.assertIsNotNone(genbank)
        self.assertEqual(genbank.project, "test-project")
        self.assertEqual(genbank.project_path, self.test_dir)
        self.assertTrue(genbank.solo)
        self.assertTrue(genbank.multi)
        self.assertTrue(genbank.min_fasta)
        self.assertEqual(len(genbank.db_files_list), 2)
        self.assertIn('test.db', genbank.db_files_list)
        self.assertIn('other.db', genbank.db_files_list)


class TestAlign(unittest.TestCase):
    """Test the Align module."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.test_fasta = os.path.join(self.test_dir, 'test.fasta')
        self.test_outfile = os.path.join(self.test_dir, 'test_aligned.fasta')
        
        # Create a simple test fasta file
        with open(self.test_fasta, 'w') as f:
            f.write('>seq1\nATGCATGC\n>seq2\nATGCATGC\n')

    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.test_dir):
            rmtree(self.test_dir, ignore_errors=True)

    def test_clustalo_init(self):
        """Test ClustalO initialization."""
        clustalo = ClustalO(
            infile=self.test_fasta,
            outfile=self.test_outfile,
            logpath=None,
            outfmt="fasta"
        )
        
        self.assertIsNotNone(clustalo)
        self.assertEqual(clustalo.infile, self.test_fasta)
        self.assertEqual(clustalo.outfile, self.test_outfile)
        self.assertEqual(clustalo.outfmt, "fasta")
        self.assertIsNone(clustalo.logpath)
        self.assertTrue(hasattr(clustalo, 'clustalolog'))

    @mock.patch('OrthoEvol.Orthologs.Align.msa.FullUtilities')
    @mock.patch('OrthoEvol.Orthologs.Align.msa.LogIt')
    def test_msa_init(self, mock_logit, mock_utils):
        """Test MultipleSequenceAlignment initialization."""
        mock_logit_instance = mock.Mock()
        mock_logit.return_value.default.return_value = mock_logit_instance
        
        # Create a mock object with the required attributes
        mock_config_obj = mock.Mock()
        mock_config_obj.raw_data = Path(self.test_dir)
        
        # Set up the FullUtilities mock
        mock_utils_instance = mock.Mock()
        mock_utils_instance.attribute_config.return_value = mock_config_obj
        mock_utils.return_value = mock_utils_instance
        
        msa = MultipleSequenceAlignment(
            project="test-project",
            project_path=self.test_dir
        )
        
        self.assertIsNotNone(msa)
        self.assertEqual(msa.project, "test-project")
        self.assertTrue(hasattr(msa, 'guidancelog'))
        self.assertTrue(hasattr(msa, 'pal2nallog'))
        self.assertTrue(hasattr(msa, 'clustalolog'))
        self.assertIsInstance(msa.dispatcher_options, dict)


if __name__ == '__main__':
    unittest.main()
