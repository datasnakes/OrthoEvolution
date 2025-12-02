import unittest
from OrthoEvol.Cookies import CookBook, Oven
from pathlib import Path
import os
import shutil

class TestCookBook(unittest.TestCase):

    def test_init(self):
        """Test CookBook initialization."""
        cookbook = CookBook()
        self.assertTrue(hasattr(cookbook, 'CookieJar'))
        self.assertIsInstance(cookbook.CookieJar, Path)
        self.assertTrue(cookbook.CookieJar.exists())

class TestOven(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures."""
        self.cookbook = CookBook()
        self.oven = Oven(recipes=self.cookbook)
        self.test_dir = Path('test_cookies_output')
        # Clean up any existing test directory
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
        self.test_dir.mkdir()
        # Clean up any leftover Test-Repository from previous runs
        test_repo = Path('Test-Repository')
        if test_repo.exists():
            shutil.rmtree(test_repo)

    def tearDown(self):
        """Clean up test directories."""
        # Clean up test directory
        if self.test_dir.exists():
            try:
                shutil.rmtree(self.test_dir)
            except (OSError, PermissionError):
                # If cleanup fails, try to remove contents individually
                for item in self.test_dir.iterdir():
                    try:
                        if item.is_dir():
                            shutil.rmtree(item)
                        else:
                            item.unlink()
                    except (OSError, PermissionError):
                        pass
                # Try to remove the directory again
                try:
                    self.test_dir.rmdir()
                except (OSError, PermissionError):
                    pass
        # Clean up any Test-Repository that might have been created
        test_repo = Path('Test-Repository')
        if test_repo.exists():
            try:
                shutil.rmtree(test_repo)
            except (OSError, PermissionError):
                pass

    def test_bake_the_repo(self):
        """Test that bake_the_repo creates a repository."""
        self.oven.repo = 'test_repo'
        self.oven.bake_the_repo(cookie_jar=self.test_dir)
        expected_dir = self.test_dir / 'test_repo'
        self.assertTrue(expected_dir.exists())

    def test_bake_the_user(self):
        """Test that bake_the_user creates a user directory."""
        self.oven.user = 'test_user'
        self.oven.bake_the_user(cookie_jar=self.test_dir)
        expected_dir = self.test_dir / 'test_user'
        self.assertTrue(expected_dir.exists())

if __name__ == '__main__':
    unittest.main()
