"""This is the test suite for Tools."""

import unittest
from shutil import rmtree

from Datasnakes.Tools import LogIt, Multiprocess



class TestTools(unittest.TestCase):
    """Test the Tools module."""

    def setUp(self, project='test-project', repository='test-repo'):
        self.project = project
        self.repo = repository

    def delete_repo(self):
        rmtree(self.repo)

    def test_logit(self):
        print()

    def test_multiprocess(self):
        print()


if __name__ == '__main__':
    unittest.main()