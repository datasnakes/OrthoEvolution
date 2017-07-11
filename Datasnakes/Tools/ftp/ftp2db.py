"""This class takes advantage of python's ftplip module and makes it easier to
download desired databases or files from NCBI."""
# Modules Used
import os
import pandas as pd
from ftplib import FTP
import time
import shutil
import subprocess
import configparser
#from cursesmenu import *  # for Linux Only
#from cursesmenu.items import *  # for Linux Only

class FTP2DB(object):
    global user
    ncbi = configparser.ConfigParser()
    ncbi.read('ncbiftp.cfg')
    """Private variable initialization"""
    __NCBI_FTP = ncbi['FTPSITE']['ncbi']
    __NCBI_RSYNC = 'rsync://ftp.ncbi.nlm.nih.gov'
    __class_name = ''
    __update_dict = {}
    __server_dict = {}
    __user = user
    __path_list = []  # Uses user input in order to change the desired self.path
    __file_list = []  # A list of files in the particular directory
    __dir_list = []  # A dynamic list of current directories
    __log_files = []
    __download_flag = False
    __ftp_path_ui = True
    __table = None

    def __init__(self, email='', path='/', ftp_update=False, db_update=False):

        # Attributes based on Ftp2Db class parameters
        self.__class_name = str(self.__class__.__name__)
        self.email = email
        self.path = path
        # Only use if there is a proper log file to be handled
        self.ftp_update_flag = ftp_update
        # Only use if there is a proper log file to be handled
        self.db_update_flag = db_update
        self.__table = [
            'blast',
            'genbank',
            'fasta'
        ]
        # Initialize a list of log files that already exist
        for item in os.listdir(where.LOG):
            if str(self.__class_name).lower() in str(item).lower():
                self.__log_files.append(item)

        if (self.ftp_update_flag is False) and (self.db_update_flag is False):
            # Other attributes
            self.database_flag = self.u_i("db_ui", table=self.__table)
            self.extensions = []  # A list of extensions in the chosen directory
            self.ext_choice = []  # The user chooses which extensions to use
            self.file_choice = []  # The user chooses which file types to use

            if self.path != '/':  # Begins Download inside the given path
                print('\nYou must already know where you want to go on NCBI\'s FTP site...'
                      '\nTaking steps to prepare for file transfer...\n')
                self.ftp_download(self.path)
            else:                 # Begins Download after navigating to the desired path
                self.ftp_navigate()

        self.update()  # ALL FILES ARE UPDATED HERE
        # Initializes __update_dict with class variables if ftp_update_flag is true
        # If db_update_flag is true it updates the databases

    def ftp_check(self):
        """Check to see if the FTP connection still exists.

        If it doesn't then it reconnects and returns an instance of the
        connection.
        """
        ftp = self.ftp_connect(self.__NCBI_FTP, self.email)
        ftp.voidcmd('NOOP')
        ftp.cwd(self.path)
        return ftp

    @staticmethod
    def ftp_connect(ftpsite, email):
        """Connect to the FTP server."""
        ftp = FTP(ftpsite, timeout=None)
        ftp.login(user='anonymous', passwd=email)
        return ftp

    @staticmethod
    def ftp_mkdir(path, p_list, update):
        """If the ftp or db files are intended to be updated then archive
        any old files. If this is a new download then make the appropriate
        directories. """
        if update is True:
            os.system('rsync -av --delete ')
            path = where.dir_archive(path, p_list)
        else:
            path = where.dir_make(path, p_list)
        return path

    @staticmethod
    def ftp_unzip(local_dir, downloaded_list):
        """Create a new unzipped file.
        Deletes the old compressed file that was downloaded."""
        for f in os.listdir(local_dir):
            _p = local_dir + '/' + f
            if f.endswith(".tar.gz"):
                # Unzip the database files
                os.system("do tar xvf " + _p + "; done")
                os.system("rm -r " + _p)
                print('Unzipped %s' % f)
            elif f.endswith(".tar"):
                os.system("tar xvf " + _p)  # Unzip the database files
                os.system("rm -r " + _p)
                print('Unzipped %s' % f)
            elif f.endswith(".gz"):
                os.system("gunzip " + _p)  # Unzip the database files
                print('Unzipped %s' % f)