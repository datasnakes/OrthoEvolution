"""File/folder zipping utility."""
import zipfile
import os
import shutil
from pathlib import Path


class ZipUtils:
    """The ZipUtils class allows easy compression/zipping of file folders.

    Inspired by http://stackoverflow.com/a/670635/7351746
    """

    def __init__(self, comp_filename, comp_path, description=None):
        """Initialize the input files and path.

        :param comp_filename (string):  This is the name of the compressed file that will be generated (eg 'test.zip')
        :param comp_path: This is the absolute path of the directory (or file) to be zipped.
        :returns:  A zip file that is created inside of the zip_path.  The path string is returned.
        """
        self.comp_filename = Path(comp_filename)
        self.comp_path = Path(comp_path)
        self.paren_comp_path = self.comp_path.parent
        self.abs_comp_file_path = self.paren_comp_path / self.comp_filename

        self.readme_path = self.comp_path / Path("README." + str(self.comp_filename) + ".txt")
        self.ignore_parts = self.comp_path.parent.parts

    def compress(self, compression=zipfile.ZIP_LZMA):
        """Zip a folder."""
        zip_handle = zipfile.ZipFile(self.abs_comp_file_path, 'w', compression)
        if os.path.isfile(self.comp_path):
            zip_handle.write(self.comp_path)
        else:
            print('skipped')
            with open(self.readme_path, 'w') as readme:
                ts = self._get_size(units="KB")
                readme.write("Compression of %s which contains %s amount of data." % (self.comp_path, ts))
            shutil.copy(str(self.readme_path), str(Path(self.comp_path).parent))
            self.__addfolder2zip(zip_handle, self.comp_path)
        zip_handle.close()
        return self.abs_comp_file_path

    def decompress(self, path=os.getcwd(), file=None):
        if not file:
            file = self.abs_comp_file_path
        zip_handle = zipfile.ZipFile(file=str(file))
        Path(path).mkdir()
        zip_handle.extractall(path=path)

    def _get_size(self, units="B"):
        total_size = 0
        unit_options = {
            "B": 1,
            "KB": 1024,
            "MB": 1048576,
            "GB": 1073741824,
            "TB": 1099511627776
        }
        for dirpath, dirnames, filenames in os.walk(self.comp_path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                total_size += os.path.getsize(fp)
        total_size = (total_size / unit_options[units])
        total_size = "%s %s" % (str(total_size), units)
        return total_size

    def __delete_contents(self, folder):
        pass

    def __addfolder2zip(self, zip_handle, folder):
        """Not meant to be used explicitly.  Use compress.

        :param zip_handle: An initialized zipfile.ZipFile handle.
        :param folder: A path that represents an entire folder to be zipped.
        :return: Recursively zips nested directories.
        """
        for file in os.listdir(folder):
            full_path = os.path.join(folder, file)
            rel_path = Path(full_path)
            rel_path = rel_path.relative_to(Path(self.comp_path))
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
