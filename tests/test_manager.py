"""This is the test suite for Manager."""
import unittest
from shutil import rmtree

from OrthoEvol.Manager.management import ProjectManagement


class TestManager(unittest.TestCase):
    """Test the Manager module."""

    def setUp(self, project='test-project', repository=None):
        self.project = project
        self.repo = repository

    def delete_project(self):
        rmtree(self.project)

    def test_projectmanagement(self):
        """Test the ProjectManagement class."""
        ProjectManagement(repo=self.repo, user=None,
                          project=self.project,
                          research=None,
                          research_type='comparative_genetics',
                          new_repo=False, new_user=False, new_project=True,
                          new_research=False)
        self.assertEqual(str(self.project), 'test-project')
        self.delete_project()


if __name__ == '__main__':
    unittest.main()
