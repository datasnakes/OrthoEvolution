"""Test the ftp module. This won't remain."""
from ftplib import FTP, error_perm
import os
import fnmatch


class NCBIFTPTest(object):
    # TODO write a real and better test.
    """A simple example for ftp."""
    def __init__(self, passwd, user='', site='ftp.ncbi.nlm.nih.gov/',
                 path='blast'):

        # Connect to the NCBI directory
        self.ftpsite = site
        sitepath = path
        ftp = FTP(self.ftpsite, timeout=None)

        # Login using email as password
        ftp.login(user=user, passwd=passwd)

        # Change to the desired directory
        # ftp.cwd(refseqrna)
        ftp.cwd(sitepath)
        # Use ftp.pwd() to find out the current directory

        # Get a list of all the files in the directory
        # ftp.retrlines('LIST')

        # This is a list of the file names
        filenames = ftp.nlst()
        file2download = 'blastftp*'

        # Create a for loop that downloads the files
        for filename in filenames:
            if fnmatch.fnmatch(filename, file2download):
                host_file = os.path.join(filename)
                try:
                    with open(host_file, 'wb') as local_file:
                        ftp.retrbinary('RETR %s' % filename, local_file.write)
                except error_perm:
                    pass

        ftp.quit()

NCBIFTPTest(passwd='shutchins2@umc.edu')