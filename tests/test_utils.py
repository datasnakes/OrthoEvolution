import unittest
import os
import time
import shutil
import tempfile

from OrthoEvol.utilities import FunctionRepeater, CookieUtils


class TestUtils(unittest.TestCase):
    """Test the Tools module."""

    def test_function_repeater(self):
        # Create a dummy function that increments a counter
        counter = 0
        def increment():
            nonlocal counter
            counter += 1
        
        # Create a FunctionRepeater instance that runs the dummy function every 1 second
        repeater = FunctionRepeater(1, increment)
        
        # Wait for 2 seconds
        time.sleep(5)
        
        # Stop the repeater
        repeater.stop()
        
        # Check that the dummy function was run at least twice
        assert counter >= 2

    def test_archive(self):
        # Create a temporary directory to use as the database path
        database_path = tempfile.mkdtemp()
        # Create a file in the database path
        with open(os.path.join(database_path, "test.txt"), "w") as f:
            f.write("Test content")
        # Create a temporary directory to use as the archive path
        archive_path = tempfile.mkdtemp()
        
        # Create a CookieUtils instance
        utils = CookieUtils()
        
        # Run the archive method with the "Full" option
        utils.archive(database_path, archive_path, "Full")
        
        # Check that the file has been archived in a .tar.xz file
        archived_files = os.listdir(archive_path)
        assert len(archived_files) == 1
        assert archived_files[0].endswith(".tar.xz")
        
        # Clean up
        shutil.rmtree(database_path)
        shutil.rmtree(archive_path)


if __name__ == '__main__':
    unittest.main()
