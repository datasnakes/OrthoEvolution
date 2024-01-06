import unittest
from OrthoEvol.Cookies import CookBook, Oven
from pathlib import Path
import os

class TestCookBook(unittest.TestCase):

    def test_init(self):
        cookbook = CookBook()
        self.assertTrue(hasattr(cookbook, 'CookieJar'))
        self.assertTrue(isinstance(cookbook.CookieJar, Path))
        # Test other attributes similarly

    def test_new_recipes(self):
        new_recipe_path = Path('path/to/new/recipe')
        cookbook = CookBook(new_recipe='new_recipe_path')
        self.assertEqual(cookbook.new_recipe, new_recipe_path)

class TestOven(unittest.TestCase):

    def setUp(self):
        self.cookbook = CookBook()
        self.oven = Oven(recipes=self.cookbook)
        self.test_dir = Path('test_directory')
        if not self.test_dir.exists():
            os.makedirs(self.test_dir)

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
        self.assertTrue((self.test_dir / repo_name).exists())

    def test_bake_the_user(self):
        user_name = 'test_user'
        self.oven.user = user_name
        self.oven.bake_the_user(cookie_jar=self.test_dir)
        self.assertTrue((self.test_dir / user_name).exists())

    # Similar tests for other methods like bake_the_project, bake_the_db_repo, etc.

if __name__ == '__main__':
    unittest.main()
