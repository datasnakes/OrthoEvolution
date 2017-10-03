"""This is the test suite for Tools."""

import unittest
import os

from Datasnakes.Tools import LogIt, Multiprocess


class TestTools(unittest.TestCase):
    """Test the Tools module."""

    def setUp(self, logfile='test.log'):
        self.logfile = logfile

    def test_logit(self):
        logit = LogIt()
        test = logit.default('testlog', self.logfile)
        self.assertEqual(str(test.name), 'TESTLOG')
        self.assertTrue(os.path.isfile(self.logfile))
        logit.shutdown()
        logit.deletelog(self.logfile)

    def test_multiprocess(self):
        print()


if __name__ == '__main__':
    unittest.main()
