"""s2s sets up sending files to servers via public SSH keys."""
import os
from pathlib import Path
import zipfile
import subprocess

# TODO-ROB:  Create Command 2 Server class or find something on GitHub similar
class S2S(object):
    """S2S (Send 2 Server) is designed for use with a public ssh key."""
    # TIP Create a public key to use this class. It's easy!
    # TIP Go here if Linux >>> http://tinyurl.com/pccz3pj

    def __init__(self, username=None, server_address=None, dest_path=None,
                 remote=True, comp_filename='', zip_path=None, compressed=True, auto=False):
        """
        :param username (string):  Remote server username.
        :param server_address (string):  Remote server address.
        :param dest_path (string):  Remote server destination.
        :param comp_filename (string):  This is the name of the compressed file that will be generated (eg 'test.zip')
        :param zip_path: This is the absolute path of the directory (or file) to be zipped.
        :param auto (bool):  Flag for automatically carrying out compression and SCP file sending.
        :param remote (bool):  Flag for determining if the copied file is on a remote or local server.  Default to True.
        """

        self.user = username
        self.address = server_address
        self.dest_path = dest_path
        self.comp_filename = comp_filename
        self.zip_path = zip_path
        self.send_cmd = "scp %s %s@%s:%s" % (self.comp_filename, self.user,
                                             self.address, self.dest_path)
        # TODO-ROB Fix this for recursive directory or just a file
        self.copy_cmd = 'cp -R %s/%s %s' % (self.zip_path,
                                            self.comp_filename, self.dest_path)
        if compressed:
            self.to_zip()
            if auto:
                if remote:
                    self.scpto(self.comp_filename)
                if not remote:
                    self.cpto(self.comp_filename)

    def scpto(self, file):
        """Send the file."""
        cmd = self.send_cmd
        status = subprocess.call([cmd], shell=True)
        if status == 0:  # Command was successful.
            print("%s file sent." % file)
        else:  # Unsuccessful. Stdout will be '1'.
            print("%s file not sent." % file)

    def cpto(self, file):
        cmd = self.copy_cmd
        status = subprocess.call([cmd], shell=True)
        if status == 0:
            print("%s file sent." % file)
        else:
            print("%s file not sent." % file)

    def to_zip(self):
        # XXX Inspired by http://tinyurl.com/y7uxn2dk
        comp_path = os.path.join(self.zip_path, self.comp_filename)
        zip_handle = zipfile.ZipFile(comp_path, 'w', zipfile.ZIP_DEFLATED)
        if os.path.isfile(self.zip_path):
            zip_handle.write(self.zip_path)
        else:
            print('skipped')
            self.__addfolder2zip(zip_handle, self.zip_path)
        zip_handle.close()
        return comp_path

    def __addfolder2zip(self, zip_handle, folder):
        """Not meant to be used explicitly.  Use to_zip."""
        # XXX Use to_zip !!!
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
