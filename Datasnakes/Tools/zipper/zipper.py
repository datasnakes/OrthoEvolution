"""File/folder zipping utility."""
import zipfile
import os
from pathlib import Path


class ZipUtils:
    """The ZipUtils class allows easy compression/zipping of file folders.

    Inspired by http://stackoverflow.com/a/670635/7351746
    """

    def __init__(self, zip_filename, zip_path):
        """Initialize the input files and path.

        :param comp_filename (string):  This is the name of the compressed file that will be generated (eg 'test.zip')
        :param zip_path: This is the absolute path of the directory (or file) to be zipped.
        :returns:  A zip file that is created inside of the zip_path.  The path string is returned.
        """
        self.zip_filename = zip_filename
        self.zip_path = zip_path
        self.ignore_parts = Path(zip_path).parent.parts

    def to_zip(self):
        """Zip a folder."""
        comp_path = os.path.join(self.zip_path, self.zip_filename)
        zip_handle = zipfile.ZipFile(comp_path, 'w', zipfile.ZIP_DEFLATED)
        if os.path.isfile(self.zip_path):
            zip_handle.write(self.zip_path)
        else:
            print('skipped')
            self.__addfolder2zip(zip_handle, self.zip_path)
        zip_handle.close()
        return comp_path

    def __addfolder2zip(self, zip_handle, folder):
        """Not meant to be used explicitly.  Use to_zip.

        :param zip_handle: An initialized zipfile.ZipFile handle.
        :param folder: A path that represents an entire folder to be zipped.
        :return: Recursively zips nested directories.
        """
        for file in os.listdir(folder):
            full_path = os.path.join(folder, file)
            rel_path = Path(full_path)
            rel_path = rel_path.relative_to(Path(self.zip_path))
            if os.path.isfile(full_path):
                if str(file) == str(self.comp_filename):
                    continue
                print('File added: ' + str(full_path))
                zip_handle.write(full_path, rel_path)
            elif os.path.isdir(full_path):
                if str(file) in self.ignore_parts:
                    continue
                print('Entering folder: ' + str(full_path))
                self.__addfolder2zip(zip_handle, full_path)
