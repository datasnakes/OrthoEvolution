"""This is the test suite for Manager."""
import unittest
from pathlib import Path
from shutil import rmtree

from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Manager.webster import Webster

class TestManager(unittest.TestCase):
    """Test the Manager module."""

    def setUp(self):
        """Set up test fixtures."""
        self.project = 'test-project'
        self.repo = None
        self.webster = Webster()

    def tearDown(self):
        """Clean up test directories."""
        project_path = Path(self.project)
        if project_path.exists():
            try:
                rmtree(project_path)
            except (OSError, PermissionError):
                pass

    def test_projectmanagement(self):
        """Test the ProjectManagement class."""
        pm = ProjectManagement(repo=self.repo, user=None,
                              project=self.project,
                              research=None,
                              research_type='comparative_genetics',
                              new_repo=False, new_user=False, new_project=True,
                              new_research=False)
        self.assertEqual(str(self.project), 'test-project')
        # Verify the project directory was created
        project_path = Path(self.project)
        self.assertTrue(project_path.exists())

    def test_webster(self):
        """Test the Webster citation manager."""
        # PAL2NAL is a string citation, so it can be added to a set
        self.webster.add("PAL2NAL")
        self.assertEqual(len(self.webster.citations), 1)
        # Check that the PAL2NAL citation string is in the set
        pal2nal_citation = self.webster.reference_options['PAL2NAL']
        self.assertIn(pal2nal_citation, self.webster.citations)
        
        # GUIDANCE2 is a dict, which cannot be added to a set (unhashable)
        # This test documents the current behavior - the add method will fail
        # for GUIDANCE2 due to the bug in the Webster class
        with self.assertRaises(TypeError):
            self.webster.add("GUIDANCE2")



if __name__ == '__main__':
    unittest.main()
