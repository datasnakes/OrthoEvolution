"""This is the test suite for Tools."""

import unittest
import os
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


if __name__ == '__main__':
    unittest.main()
