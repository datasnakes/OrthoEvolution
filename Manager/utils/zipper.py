import zipfile
import os
from pathlib import Path


class ZipUtilities:

    def __init__(self, comp_filename, zip_path):
        self.comp_filename = comp_filename
        self.zip_path = zip_path
        self.ignore_parts = Path(zip_path).parent.parts

    def to_zip(self):
        comp_path = os.path.join(self.zip_path, self.comp_filename)
        zip_file = zipfile.ZipFile(comp_path, 'w', zipfile.ZIP_DEFLATED)
        if os.path.isfile(self.zip_path):
            zip_file.write(self.zip_path)
        else:
            print('skipped')
            self.add_folder_to_zip(zip_file, self.zip_path)
        zip_file.close()

    def add_folder_to_zip(self, zip_file, folder):
        for file in os.listdir(folder):
            full_path = os.path.join(folder, file)
            rel_path = Path(full_path)
            rel_path = rel_path.relative_to(Path(self.zip_path))
            if os.path.isfile(full_path):
                if str(file) == str(self.comp_filename):
                    continue
                print('File added: ' + str(full_path))
                zip_file.write(full_path, rel_path)
            elif os.path.isdir(full_path):
                if str(file) in self.ignore_parts:
                    continue
                print('Entering folder: ' + str(full_path))
                self.add_folder_to_zip(zip_file, full_path)
