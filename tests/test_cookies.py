import unittest
import os

from OrthoEvol.Cookies import CookBook, Oven

class TestCookBookOven(unittest.TestCase):
    def test_cookbook_init(self):
        # Create a new instance of the CookBook class
        cookbook = CookBook()

        # Check that the attributes of the CookBook class are correctly initialized
        self.assertTrue(hasattr(cookbook, "CookieJar"))
        self.assertTrue(hasattr(cookbook, "repo_cookie"))
        self.assertTrue(hasattr(cookbook, "user_cookie"))
        self.assertTrue(hasattr(cookbook, "project_cookie"))
        self.assertTrue(hasattr(cookbook, "basic_project_cookie"))
        self.assertTrue(hasattr(cookbook, "research_cookie"))
        self.assertTrue(hasattr(cookbook, "app_cookie"))
        self.assertTrue(hasattr(cookbook, "db_cookie"))
        self.assertTrue(hasattr(cookbook, "website_cookie"))

    def test_oven_bake_cookies(self):
        # Create a new instance of the Oven class
        oven = Oven()

        # Set the output directory for the baked cookies
        cookie_jar = "test_cookie_jar"
        oven.cookie_jar = cookie_jar

        # Bake a cookie and check that a new directory was created in the output directory
        oven.bake_cookies(recipe="repo_cookie", repo="test_repo")
        self.assertTrue(os.path.exists(os.path.join(cookie_jar, "test_repo")))

if __name__ == '__main__':
    unittest.main()
