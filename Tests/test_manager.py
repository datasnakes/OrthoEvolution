"""This is the test suite for Manager."""
import unittest
from shutil import rmtree

from Datasnakes.Manager import ProjectManagement


class TestManager(unittest.TestCase):
    """Test the Manager module."""

    def setUp(self, project='test-project', repository='test-repo'):
        self.project = project
        self.repo = repository

    def delete_repo(self):
        rmtree(self.repo)

    def test_projectmanagement(self):
        """Test the ProjectManagement class."""
        ProjectManagement(repo=self.repo, user='test-user', project=self.project,
                          research='test-research', research_type='comparative_genetics',
                          new_repo=True, new_user=True, new_project=True,
                          new_research=True)
        self.assertEqual(str(self.repo), 'test-repo')
        self.delete_repo()
    # TODO add tests for each individual subclass


if __name__ == '__main__':
    unittest.main()
