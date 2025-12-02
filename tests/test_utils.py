import unittest
from pathlib import Path
import os
import shutil

from OrthoEvol.utilities import CookieUtils, FunctionRepeater


class TestCookieUtils(unittest.TestCase):

    def setUp(self):
        self.utils = CookieUtils()
        self.test_dir = Path('test_dir')
        self.test_dir.mkdir(exist_ok=True)
        self.archive_path = Path('archive_dir')
        self.archive_path.mkdir(exist_ok=True)

    def tearDown(self):
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
        if self.archive_path.exists():
            shutil.rmtree(self.archive_path)

    def test_archive(self):
        # Mocking file and directory creation for the test
        test_file = self.test_dir / 'test.txt'
        with open(test_file, 'w') as f:
            f.write('test')
        self.assertTrue(test_file.exists())

        # Test archive functionality
        archive_list = self.utils.archive(database_path=self.test_dir, archive_path=self.archive_path, option='Full')
        self.assertIsInstance(archive_list, list)
        self.assertTrue(any(self.archive_path in Path(a) for a in archive_list))

    def test_get_size(self):
        test_file = self.test_dir / 'test.txt'
        with open(test_file, 'w') as f:
            f.write('test')
        size = self.utils.get_size(start_path=str(test_file))
        self.assertIsInstance(size, str)

class TestFunctionRepeater(unittest.TestCase):

    def setUp(self):
        self.mock_function = unittest.mock.Mock()
        self.repeater = FunctionRepeater(interval=1, function=self.mock_function)

    def tearDown(self):
        self.repeater.stop()

    def test_repeater_start_stop(self):
        self.assertTrue(self.repeater.is_running)
        self.repeater.stop()
        self.assertFalse(self.repeater.is_running)

if __name__ == '__main__':
    unittest.main()
