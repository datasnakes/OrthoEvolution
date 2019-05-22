"""This is the test suite for Tools."""

import unittest
import os

from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Tools.parallel import Multiprocess
from OrthoEvol.Tools.ftp import NcbiFTPClient


class TestTools(unittest.TestCase):
    """Test the Tools module."""

    def setUp(self, logfile='test.log', filename='README.ftp'):
        self.logfile = logfile
        self.filename = filename

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

    def test_ncbiftpclient(self):
        """Test the NcbiFTPClient."""
        ncbiftp = NcbiFTPClient(self.assertTrue('someone@gmail.com'))
        self.assertIsNotNone(ncbiftp.ftp.welcome)
        ncbiftp.close_connection()


if __name__ == '__main__':
    unittest.main()
