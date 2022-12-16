"""This is the test suite for Tools."""

import unittest
import os
from pkg_resources import resource_filename

from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Tools.parallel import Multiprocess
from OrthoEvol.Tools.ftp import NcbiFTPClient
from OrthoEvol.Tools.mygene import MyGene
from OrthoEvol.Manager.config import test


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


if __name__ == '__main__':
    unittest.main()
