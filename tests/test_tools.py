"""This is the test suite for Tools."""

import unittest
from unittest import mock
import os
import logging
import tempfile
from pkg_resources import resource_filename

from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Tools.parallel import Multiprocess
from OrthoEvol.Tools.ftp import NcbiFTPClient
from OrthoEvol.Tools.mygene import MyGene
from OrthoEvol.Manager.config import test
from OrthoEvol.Tools.pybasher import PyBasher


class TestTools(unittest.TestCase):
    """Test the Tools module."""

    def setUp(self, logfile='test.log'):
        self.logfile = logfile
        self.infile = resource_filename(test.__name__, 'test_blast.csv')
        self.outfile = 'out_mygene.csv'

    def test_logit(self):
        """Test the LogIt class."""
        logit = LogIt()
        test = logit.default(logname='testlog', logfile=self.logfile)
        self.assertEqual(str(test.name), 'TESTLOG')
        self.assertTrue(os.path.isfile(self.logfile))
        logit.shutdown()
        logit.deletelog(self.logfile)

    def test_multiprocess(self):
        """Test the Multiprocess class."""
        mp = Multiprocess()
        self.assertIsNotNone(mp.cpus)
        self.assertEqual(mp.num_procs, mp.cpus - 1)

    def test_ncbiftpclient(self):
        """Test the NcbiFTPClient class."""
        ncbiftp = NcbiFTPClient(self.assertTrue('someone@gmail.com'))
        self.assertIsNotNone(ncbiftp.ftp.welcome)
        ncbiftp.close_connection()

    def test_mygene(self):
        """Test the MyGene class."""
        mg = MyGene(infile=self.infile, outfile=self.outfile)
        self.assertEqual(mg.species, 'human')
        self.assertIsNotNone(mg.accessions_list)
        self.assertIsInstance(mg.accessions_list, list)
        mg.query_mygene()
        os.remove(self.outfile)

    def test_mv(self):
        # Create a file to move
        with open("test.txt", "w") as f:
            f.write("Test content")
        
        # Create a PyBasher instance
        pybasher = PyBasher()
        
        # Move the file
        pybasher.mv("test.txt", "moved.txt")
        
        # Check that the file was moved
        assert not os.path.exists("test.txt")
        assert os.path.exists("moved.txt")
        
        # Check that the moved file has the correct contents
        with open("moved.txt", "r") as f:
            moved_content = f.read()
        assert moved_content == "Test content"
        
        # Clean up
        os.remove("moved.txt")

    def test_cp(self):
        # Create a file to copy
        with open("test.txt", "w") as f:
            f.write("Test content")
        
        # Create a PyBasher instance
        pybasher = PyBasher()
        
        # Copy the file
        pybasher.cp("test.txt", "copy.txt")
        
        # Check that the copy was made
        assert os.path.exists("copy.txt")
        
        # Check that the copy has the same contents as the original
        with open("test.txt", "r") as f:
            original_content = f.read()
        with open("copy.txt", "r") as f:
            copy_content = f.read()
        assert original_content == copy_content
        
        # Clean up
        os.remove("test.txt")
        os.remove("copy.txt")

    def test_logit_custom(self):
        """Test LogIt custom method."""
        logit = LogIt()
        custom_log = logit.custom(
            logname='customlog',
            logfile=None,
            level=logging.INFO,
            fmt='default'
        )
        self.assertIsNotNone(custom_log)
        self.assertEqual(custom_log.name, 'customlog')
        logit.shutdown()

    def test_pybasher_mkdir_touch(self):
        """Test PyBasher mkdir and touch methods."""
        test_dir = tempfile.mkdtemp()
        try:
            pybasher = PyBasher()
            test_file = os.path.join(test_dir, 'test_file.txt')
            test_subdir = os.path.join(test_dir, 'test_subdir')
            
            # Test mkdir
            pybasher.mkdir(test_subdir)
            self.assertTrue(os.path.isdir(test_subdir))
            
            # Test touch
            pybasher.touch(test_file)
            self.assertTrue(os.path.isfile(test_file))
        finally:
            import shutil
            shutil.rmtree(test_dir, ignore_errors=True)

    def test_pybasher_pwd(self):
        """Test PyBasher pwd method."""
        pybasher = PyBasher()
        pybasher.pwd()
        # pwd returns stdout, verify it's a string
        result = str(pybasher)
        self.assertIsInstance(result, str)
        self.assertTrue(len(result) > 0)

    def test_ncbiftp_pathformat(self):
        """Test NcbiFTPClient _pathformat classmethod."""
        # Test valid path format
        valid_path = '/blast/db/'
        try:
            NcbiFTPClient._pathformat(valid_path)
        except ValueError:
            self.fail("_pathformat raised ValueError for valid path")
        
        # Test invalid path format
        invalid_path = '/blast/db'  # Missing trailing slash
        with self.assertRaises(ValueError):
            NcbiFTPClient._pathformat(invalid_path)

    @mock.patch('OrthoEvol.Tools.parallel.multiprocess.Pool')
    def test_multiprocess_map_to_function(self, mock_pool):
        """Test Multiprocess map_to_function method."""
        mp = Multiprocess()
        
        def test_func(x):
            return x * 2
        
        test_list = [1, 2, 3]
        mock_pool_instance = mock.Mock()
        mock_pool.return_value.__enter__.return_value = mock_pool_instance
        
        mp.map_to_function(test_func, test_list, procs=2)
        
        # Verify pool was created and map was called
        mock_pool.assert_called_once_with(processes=2)
        mock_pool_instance.map.assert_called_once_with(test_func, test_list)


if __name__ == '__main__':
    unittest.main()
