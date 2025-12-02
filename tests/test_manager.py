"""This is the test suite for Manager."""
import unittest
from unittest import mock
from pathlib import Path
from shutil import rmtree
import tempfile

from OrthoEvol.Manager.management import (
    Management, RepoManagement, UserManagement,
    WebsiteManagement, ProjectManagement
)
from OrthoEvol.Manager.webster import Webster
from OrthoEvol.Manager.database_management import BaseDatabaseManagement

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

    def test_webster_show(self):
        """Test Webster show method."""
        self.webster.add("PAL2NAL")
        self.webster.add("CLUSTALO")
        # show() prints to stdout, so we just verify it doesn't raise an error
        try:
            self.webster.show()
        except Exception as e:
            self.fail(f"show() raised {e} unexpectedly")

    def test_webster_init(self):
        """Test Webster initialization."""
        webster = Webster()
        self.assertIsInstance(webster.citations, set)
        self.assertEqual(len(webster.citations), 0)
        self.assertIn('PAL2NAL', webster.reference_options)
        self.assertIn('Full', webster.archive_options)


class TestManagement(unittest.TestCase):
    """Test the Management base class."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Clean up test directories."""
        if self.test_dir.exists():
            try:
                rmtree(self.test_dir)
            except (OSError, PermissionError):
                pass

    def test_management_init(self):
        """Test Management initialization."""
        mgmt = Management(repo=None, home=str(self.test_dir))
        self.assertIsNotNone(mgmt.file_home)
        self.assertIsNotNone(mgmt.Kitchen)
        self.assertIsNotNone(mgmt.Manager)
        self.assertIsNotNone(mgmt.Orthologs)
        self.assertIsNotNone(mgmt.Tools)

    def test_management_with_repo(self):
        """Test Management with repository."""
        mgmt = Management(repo='test-repo', home=str(self.test_dir))
        self.assertEqual(mgmt.repo, 'test-repo')
        self.assertIsNotNone(mgmt.repo_path)


class TestRepoManagement(unittest.TestCase):
    """Test the RepoManagement class."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.repo = 'test-repo'

    def tearDown(self):
        """Clean up test directories."""
        repo_path = self.test_dir / self.repo
        if repo_path.exists():
            try:
                rmtree(repo_path)
            except (OSError, PermissionError):
                pass
        if self.test_dir.exists():
            try:
                rmtree(self.test_dir)
            except (OSError, PermissionError):
                pass

    def test_repo_management_init(self):
        """Test RepoManagement initialization."""
        repo_mgmt = RepoManagement(
            repo=self.repo,
            user=None,
            home=str(self.test_dir),
            new_repo=False,
            new_user=False
        )
        self.assertEqual(repo_mgmt.repo, self.repo)
        self.assertIsNotNone(repo_mgmt.repo_path)
        self.assertIsNotNone(repo_mgmt.docs)
        self.assertIsNotNone(repo_mgmt.users)

    def test_repo_management_with_user(self):
        """Test RepoManagement with user."""
        repo_mgmt = RepoManagement(
            repo=self.repo,
            user='test-user',
            home=str(self.test_dir),
            new_repo=False,
            new_user=False
        )
        self.assertEqual(repo_mgmt.user, 'test-user')
        self.assertIsNotNone(repo_mgmt.user_path)


class TestUserManagement(unittest.TestCase):
    """Test the UserManagement class."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.repo = 'test-repo'
        self.user = 'test-user'

    def tearDown(self):
        """Clean up test directories."""
        repo_path = self.test_dir / self.repo
        if repo_path.exists():
            try:
                rmtree(repo_path)
            except (OSError, PermissionError):
                pass
        if self.test_dir.exists():
            try:
                rmtree(self.test_dir)
            except (OSError, PermissionError):
                pass

    def test_user_management_init(self):
        """Test UserManagement initialization."""
        # UserManagement requires project to be set (even if None) for Oven initialization
        # So we pass an empty string or a dummy project name
        user_mgmt = UserManagement(
            repo=self.repo,
            user=self.user,
            project='test-project',
            home=str(self.test_dir),
            new_user=False,
            new_project=False
        )
        self.assertEqual(user_mgmt.user, self.user)
        self.assertEqual(user_mgmt.repo, self.repo)
        self.assertEqual(user_mgmt.project, 'test-project')
        self.assertIsNotNone(user_mgmt.user_path)
        self.assertIsNotNone(user_mgmt.user_db)
        self.assertIsNotNone(user_mgmt.projects)

    def test_user_management_with_project(self):
        """Test UserManagement with project."""
        user_mgmt = UserManagement(
            repo=self.repo,
            user=self.user,
            project='test-project',
            home=str(self.test_dir),
            new_user=False,
            new_project=False
        )
        self.assertEqual(user_mgmt.project, 'test-project')
        self.assertIsNotNone(user_mgmt.project_path)


class TestWebsiteManagement(unittest.TestCase):
    """Test the WebsiteManagement class."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.repo = 'test-repo'
        self.website = 'test-website'

    def tearDown(self):
        """Clean up test directories."""
        repo_path = self.test_dir / self.repo
        if repo_path.exists():
            try:
                rmtree(repo_path)
            except (OSError, PermissionError):
                pass
        if self.test_dir.exists():
            try:
                rmtree(self.test_dir)
            except (OSError, PermissionError):
                pass

    def test_website_management_init(self):
        """Test WebsiteManagement initialization."""
        web_mgmt = WebsiteManagement(
            repo=self.repo,
            website=self.website,
            home=str(self.test_dir),
            new_website=False
        )
        self.assertEqual(web_mgmt.website, self.website)
        self.assertEqual(web_mgmt.web_host, '0.0.0.0')
        self.assertEqual(web_mgmt.web_port, '5252')
        self.assertIsNotNone(web_mgmt.website_path)

    def test_website_management_custom_host_port(self):
        """Test WebsiteManagement with custom host and port."""
        web_mgmt = WebsiteManagement(
            repo=self.repo,
            website=self.website,
            host='127.0.0.1',
            port='8080',
            home=str(self.test_dir),
            new_website=False
        )
        self.assertEqual(web_mgmt.web_host, '127.0.0.1')
        self.assertEqual(web_mgmt.web_port, '8080')


class TestProjectManagementExtended(unittest.TestCase):
    """Extended tests for ProjectManagement."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.project = 'test-project-extended'
        self.repo = None

    def tearDown(self):
        """Clean up test directories."""
        project_path = Path(self.project)
        if project_path.exists():
            try:
                rmtree(project_path)
            except (OSError, PermissionError):
                pass
        if self.test_dir.exists():
            try:
                rmtree(self.test_dir)
            except (OSError, PermissionError):
                pass

    def test_project_management_with_research(self):
        """Test ProjectManagement with research."""
        pm = ProjectManagement(
            repo=self.repo,
            user=None,
            project=self.project,
            research='test-research',
            research_type='comparative_genetics',
            new_project=True,
            new_research=False
        )
        self.assertEqual(pm.research, 'test-research')
        self.assertEqual(pm.research_type, 'comparative_genetics')
        self.assertIsNotNone(pm.research_path)

    def test_project_management_with_app(self):
        """Test ProjectManagement with app."""
        pm = ProjectManagement(
            repo=self.repo,
            user=None,
            project=self.project,
            research='test-research',
            research_type='comparative_genetics',
            app='test-app',
            new_project=True,
            new_research=False,
            new_app=False
        )
        self.assertEqual(pm.app, 'test-app')
        self.assertIsNotNone(pm.app_path)

    def test_project_management_paths(self):
        """Test ProjectManagement path attributes."""
        pm = ProjectManagement(
            repo=self.repo,
            user=None,
            project=self.project,
            research=None,
            research_type='comparative_genetics',
            new_project=True,
            new_research=False
        )
        self.assertIsNotNone(pm.project_path)
        self.assertIsNotNone(pm.project_index)
        self.assertIsNotNone(pm.data)
        self.assertIsNotNone(pm.raw_data)
        self.assertIsNotNone(pm.project_web)


class TestBaseDatabaseManagement(unittest.TestCase):
    """Test the BaseDatabaseManagement class."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.email = 'test@example.com'
        self.driver = 'sqlite'

    def tearDown(self):
        """Clean up test directories."""
        if self.test_dir.exists():
            try:
                rmtree(self.test_dir)
            except (OSError, PermissionError):
                pass

    @mock.patch('OrthoEvol.Manager.database_management.NcbiFTPClient')
    def test_base_database_management_init(self, mock_ftp):
        """Test BaseDatabaseManagement initialization."""
        mock_ftp_instance = mock.Mock()
        mock_ftp.return_value = mock_ftp_instance
        db_mgmt = BaseDatabaseManagement(
            email=self.email,
            driver=self.driver,
            project='test-project',
            project_path=str(self.test_dir),
            proj_mana=None,
            ftp_flag=True
        )
        self.assertEqual(db_mgmt.email, self.email)
        self.assertEqual(db_mgmt.driver, self.driver)
        self.assertEqual(db_mgmt.project, 'test-project')
        self.assertIsNotNone(db_mgmt.database_path)

    def test_base_database_management_no_ftp(self):
        """Test BaseDatabaseManagement without FTP."""
        db_mgmt = BaseDatabaseManagement(
            email=self.email,
            driver=self.driver,
            project='test-project',
            project_path=str(self.test_dir),
            proj_mana=None,
            ftp_flag=False
        )
        self.assertFalse(db_mgmt.ftp_flag)
        self.assertFalse(hasattr(db_mgmt, 'ncbiftp'))



if __name__ == '__main__':
    unittest.main()
