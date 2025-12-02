import unittest
from unittest.mock import patch, MagicMock
from OrthoEvol.Cookies import CookBook, Oven
from pathlib import Path
import os

class TestCookBook(unittest.TestCase):

    def test_init(self):
        cookbook = CookBook()
        self.assertTrue(hasattr(cookbook, 'CookieJar'))
        self.assertIsInstance(cookbook.CookieJar, Path)

    @patch('builtins.open', new_callable=MagicMock)
    def test_new_recipes(self, mock_open):
        mock_open.return_value.__enter__.return_value = MagicMock()
        new_recipe_path = Path('path/to/new/recipe')
        cookbook = CookBook(new_recipe='new_recipe_path')
        self.assertTrue(hasattr(cookbook, 'new_recipe'))
        self.assertEqual(getattr(cookbook, 'new_recipe'), new_recipe_path)

class TestOven(unittest.TestCase):

    def setUp(self):
        self.cookbook = CookBook()
        self.oven = Oven(recipes=self.cookbook)
        self.test_dir = Path('test_directory')
        self.test_dir.mkdir(exist_ok=True)

    def tearDown(self):
        if self.test_dir.exists():
            os.rmdir(self.test_dir)

    def test_init(self):
        self.assertEqual(self.oven.Recipes, self.cookbook)
        self.assertEqual(self.oven.cookie_jar, os.getcwd())

    def test_bake_the_repo(self):
        repo_name = 'test_repo'
        self.oven.repo = repo_name
        self.oven.bake_the_repo(cookie_jar=self.test_dir)
        expected_dir = self.test_dir / repo_name
        self.assertTrue(expected_dir.exists())

    def test_bake_the_user(self):
        user_name = 'test_user'
        self.oven.user = user_name
        self.oven.bake_the_user(cookie_jar=self.test_dir)
        expected_dir = self.test_dir / user_name
        self.assertTrue(expected_dir.exists())

    # Similar tests for other methods like bake_the_project, etc.

if __name__ == '__main__':
    unittest.main()
