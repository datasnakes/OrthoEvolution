<<<<<<< HEAD
""""""
=======
# -*- coding: utf-8 -*-
"""
Orthologs-Project
Ftp2Db updated on 11/18/2016 at 3:16 PM

Input:  The class takes an email and an NCBI FTP site path
(ex..  refseq/release/vertebrate_mammalian) as optional attributes.
User input is also needed.

Output:  The user gets to choose which files to download.  A cache is also
created in order update the same files at a later date.

Description:  This class is meant to act as an interface for downloading files
from NCBI.  Calling the class begins a series questions posed to the user
to navigate the NCBI FTP site and to begin downloading files from the site.

@author: Rob Gilmore
"""
>>>>>>> b4e6bb4ddfa7bc087f1fae5a8844594a6a6198c4
# Modules Used
import os
import pandas as pd
from ftplib import FTP
import time
import shutil
import subprocess
import configparser
<<<<<<< HEAD
=======
#from cursesmenu import *  # for Linux Only
#from cursesmenu.items import *  # for Linux Only

from Datasnakes.Orthologs.CompGenetics.comp_gen import CompGenAnalysis
from Datasnakes.Manager.utils.mana import ProjMana
#class FTP2DB(object):
# Set up directories and project
home = os.getcwd()
project = "GPCR-Orthologs-Project"
user = "rgilmore"
#where = ProjMana(home, project)
# Use lister() class here so that we can easily access our Master RNA
# Accession File
# Always make sure this file name is correct
#what = CompGenAnalysis('Master_RNA_Accession_File.csv')
>>>>>>> b4e6bb4ddfa7bc087f1fae5a8844594a6a6198c4


class Ftp2Db(object):
    global user
    ncbi = configparser.ConfigParser()
    ncbi.read('ncbiftp.cfg')
    """Private variable initialization"""
    # __NCBI_FTP = ncbi['FTPSITE']['ncbi']
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

    def u_i(self, ui_value, table=None, prefix=None, suffix=None):
        """Handles user input for other instances. The ui_value determines the'
        proper output. """
        ext_choice = []
    # ftp_path_ui
    # download_ui
        if ui_value == 'ftp_path_ui':
            self.__ftp_path_ui = True
            while self.__ftp_path_ui is True:
                answer = (input('What do you want to do?\nTo move forward in the directories type the name of'
                                ' the desired directory. (BE EXACT!!!)\nTo move backwards in the directories Press '
                                '\'ENTER\'.\nTo download files type \'DOWNLOAD\''))
                # // TODO-ROB use curses module to use arrow keys for scrolling through a list of directory names
                if answer.lower() == 'download':  # Preparing to download some files (exiting u_i)
                    print('Preparing for download.')
                    self.__download_flag = True
                    self.__ftp_path_ui = False
                elif answer.lower() == '':  # Moving backwards in a directory (exiting u_i)
                    print('moving backwards')
                    if len(self.__path_list) > 0:
                        self.__path_list.pop(-1)
                    else:
                        continue
                    self.__ftp_path_ui = False
                elif answer.lower() in map(str.lower, self.__dir_list):  # Moving forward in a directory (exiting u_i)
                    # //TODO-ROB Allow for the user to put in a path with '/' as well
                    print('moving forward')
                    self.__ftp_path_ui = False
                    self.__path_list.append(answer)
                # The user didn't properly type something or something is wrong
                # (restating the u_i)
                else:
                    print('\nBe careful with your input.\n\nTry Again.\n')
                    time.sleep(3)
                    print(table.to_string(index=False, header=False) + '\n\n')
            self.__ftp_path_ui = True

    # ext_ui
        if ui_value == 'ext_ui':
            ext_dict = {}
            file_type = []
            print('\n')
            if prefix is not None:
                ext_dict['prefix'] = []
                for ft in prefix:
                    if ft == '':
                        continue
                    print('File type: ' + ft)
                    ext_dict['prefix'].append(ft)
                print('\n')
            if suffix is not None:
                ext_dict['suffix'] = []
                for extension in suffix:
                    if extension == '':
                        continue
                    print('Extensions: ' + extension)
                    ext_dict['suffix'].append(extension)

            for key, value in ext_dict.items():
                count = 0
                if str(key) == 'prefix':
                    for f_type in value:
                        if f_type is None:
                            continue
                        count += 1
                        if count == 1:
                            answer = input('\n\nPress \'ENTER\' to keep the following file types in the download que.'
                                           '\nType \'REMOVE\' to remove this type of file from the download que.\n\n %s'
                                           % f_type)
                        else:
                            answer = input('%s' % f_type)
                        if 'remove' in answer.lower():
                            continue
                        file_type.append(f_type)
                if str(key) == 'suffix':
                    for f_type in value:

                        count += 1
                        if count == 1:
                            answer = input('\n\nPress \'ENTER\' to keep the following file types in the download que.'
                                           '\nType \'REMOVE\' to remove this type of file from the download que.\n\n %s'
                                           % f_type)
                        else:
                            answer = input('%s' % f_type)
                        if 'remove' in answer.lower():
                            continue
                        ext_choice.append(f_type)
            return file_type, ext_choice
        if ui_value == 'db_ui':
            while_flag = False
            while while_flag is False:
                answer = input('What type of database file are you downloading from the NCBI-FTP site?\n'
                               'blast, genbank, or fasta?\n'
                               'If none of these apply \'Press Enter\'')

                if answer.lower() in table:
                    while_flag = True
                    return answer
                else:
                    print('Choice not available...\nType carefully next time...')

    def ftp_navigate(self):
        """Continues navigation until the user decides to download files."""
        self.__path_list = []
        if self.__download_flag is False:
            # Check the connection
            ftp = self.ftp_check()
            print(ftp.getwelcome())
            input('Press Enter')
        while self.__download_flag is False:
            print('download_flag is %s' % self.__download_flag)
            print('ui is %s' % self.__ftp_path_ui)
            self.ftp_path()
            self.path = '/'
        self.__ftp_path_ui = False
        self.ftp_path()
        self.ftp_download(self.path)

    def ftp_path(self):
        """Using user input, determine what directories to navigate on the
        NCBI FTP server."""

        # (re-initialize self.path)
        # Check the connection and change the directory
        if self.path == '/':                # *
            for item in self.__path_list:   # *
                self.path += item + '/'     # *
        ftp = self.ftp_check()              # *
        if len(self.__path_list) != 0:      # *
            ftp.cwd(self.path)              # *

        ui_value = 'ftp_path_ui'
        self.__dir_list = []
        nest_list = []
        temp_list = []

        # Display of the current items in the directory
        count = 0
        if self.__ftp_path_ui is True:
            print('\nHere is a list of what\'s in the current directory:\n')
            self.__dir_list = ftp.nlst()
            self.__dir_list.sort()
            dir_count = len(self.__dir_list)

            for item in self.__dir_list:
                temp_list.append(item + '  ')
                if len(temp_list) % 7 == 0:
                    count += 7
                    nest_list.append(temp_list)
                    temp_list = []

            if (dir_count - count) != 0:
                nest_list.append(temp_list)
                # temp_list = []
            table = pd.DataFrame(nest_list)
            print(table.to_string(index=False, header=False) + '\n\n')

            self.u_i(ui_value, table)  # User Input

    def ftp_download(self, path):
        """Initiates the downloading process which involves the user
        selecting the proper files."""

        self.__file_list = []
        # Check the connection and make a file list
        ftp = self.ftp_check()  # *
        ftp.retrlines('NLST', self.__file_list.append)

        # ui_value = 'download_ui'
        download_list = []  # Keeps track of what the user wants to download
        downloaded_list = []  # Keeps track of download-ED files
        old_archive = []  # For Logging

        # Creates a helpful naming convention for the log file based on the
        # NCBI-FTP path
        log_name = str(self.path)
        log_name = log_name[1:-1].replace('/', '-')
        self.__update_dict['log_name'] = log_name
        log_file_name = where.LOG + ('/%s_' %
                                     self.__class_name) + log_name + '.log'
        if self.ftp_update_flag is False:
            f_handle = 'w'  # Overwrites the log file when we are not updating
            self.file_choice, self.ext_choice = self.ftp_extension(
                path)  # Get the users extension choice
            local_dir, t = self.ftp_mkdir(
                where.NCBI_DATA, self.__path_list, self.ftp_update_flag)  # Handles directories

        else:
            f_handle = 'a'  # For FTP updates this makes the log file append the log file
            self.__path_list = where.path_list_make(self.path, where.NCBI_DATA)
            local_dir, t = self.ftp_mkdir(
                where.NCBI_DATA, self.__path_list, self.ftp_update_flag)  # Directory Handling
            print('\n\nThis is an update, so the previous selections will be downloaded\n'
                  'and the old files will be archived.\n\n')

            # Save the old archive line to a list
            with open(log_file_name, 'r') as log_file:
                for line in log_file:
                    if line.startswith('<<'):
                        old_archive.append(line)

            # Update the logging file and add the old/new archive lines
            with open(log_file_name, 'w') as log_file:
                for item in old_archive:
                    log_file.write(item)
                if t is not None:  # Logs the files that were archived during the current FTP session
                    log_file.write(
                        '\n<< Archived ,%s, stamped as ,%s.\n' %
                        (local_dir, t))
                    log_file.write(time.strftime("<< %m/%d/%Y %I:%M:%S\n"))

        # Make a list of files to download
        for file in self.__file_list:
            for file_choice in self.file_choice:
                if file_choice not in file:
                    continue
                else:
                    for ext_choice in self.ext_choice:
                        if ext_choice not in file:
                            continue
                        download_list.append(file)
        # TODO-ROB Make a logging function and look up the built in logging method
        # TODO-ROB Add Code that sends Rob and Shaurita a login message about the current update
        # TODO-SDH Add Code that sends Rob and Shaurita a login message about the current update

        with open(log_file_name, f_handle) as log_file:

            """
            ######################################################
            # ('> ' means that this line is code to initiate)    #
            # ('< ' means that this line is documentation)       #
            # ('# ' means that this line is a file name)         #
            # ('<< ' means tha this line is archive information) #
            ######################################################
            """

            for key, value in self.__dict__.items():
                if isinstance(
                        value, bool):  # Skips the update flags so the configuration doesn't break
                    continue
                # Log script to be used ass class variables later (a.k.a.
                # "self.key = value")
                if '__' in str(key):
                    continue
                if isinstance(value, str):
                    log_file.write(
                        '> update_dict[\'{0}\'] = \'{1}\'\n'.format(
                            key, value))
                    continue
                log_file.write(
                    '> update_dict[\'{0}\'] = {1}\n'.format(
                        key, value))  # Logs the objects in
                # order to make the updates work
            t = time.strftime("%m/%d/%Y")
            log_file.write(
                '\n<%s> Below are the files that have been downloaded and Unzipped to %s: \n' % (t, local_dir))

            for file in download_list:  # Starts the download process
                file = str(file)
                p = local_dir + '/' + file
                # Check the connection
                # *************************
                ftp = self.ftp_check()  # *
                # *************************
                try:  # If the file doesn't download all the way from keyboard interrupt then delete it
                    with open(p, 'wb') as local_file:
                        ftp.retrbinary('RETR %s' % file, local_file.write)
                except KeyboardInterrupt:
                    os.remove(local_dir + '/' + file)
                    # unzip the files that were downloaded
                    self.ftp_unzip(local_dir, downloaded_list)
                    raise KeyboardInterrupt
                t = time.strftime("%I:%M:%S")
                log_file.write(('#%s# ' % t) + file[:-3] + '\n')
                # Keeping track of the downloaded files
                downloaded_list.append(file)
                print(file + ' Created')

        self.ftp_unzip(local_dir, downloaded_list)

        if self.db_update_flag is False:
            self.db_upload(('%s_' % self.__class_name) + log_name + '.log')
        # //TODO-ROB come up with a better naming convention thats good for UPDATE module and Mana()
        # //TODO-ROB append a file with different paths for the databases for use in Mana()
        # //TODO-ROB make a logging function

        if self.ftp_update_flag is False:
            # //TODO-ROB Test this
            answer = input('All files in the chose directory have been downloaded without interruption.\n'
                           'If you want to download more files type \'DOWNLOAD\'')
            if answer.lower() == 'download':
                self.__init__(self.email, path='/')

    def ftp_extension(self, path):
        """Creates a list of extensions in the path and then presents them to
        the user for selection."""
        # Initialize globals and locals
        ui_value = 'ext_ui'
        suffix = []
        prefix = []
        temp_dict = {}
        self.__file_list = []
        count1 = 0

        # Check the connection
        ftp = self.ftp_check()

        ftp.cwd('/' + path)
        ftp.retrlines('NLST', self.__file_list.append)
        for item2 in self.__file_list:
            count1 += 1
            prefix_list = []
            suffix_list = []
            ext_list = item2.split(sep='.')
            count2 = 0
            for item1 in ext_list:
                count2 += 1
                if count2 == 1:
                    prefix_list = []
                    prefix_list.append(item1)
                if str(item1).isdigit() is True:
                    try:
                        if str(ext_list[count2]).isdigit() is True:
                            continue
                        elif str(ext_list[count2]).isdigit() is False:
                            suffix_list = ext_list[count2:(len(ext_list))]
                    except IndexError:
                        continue

            ext1 = '.'.join(suffix_list)
            ext2 = prefix_list[0]
            if count1 == 1:
                prefix.append(ext2)
                self.__server_dict[str(ext2)] = []
            if ext1 not in suffix:
                suffix.append(ext1)
            if ext2 not in prefix:
                prefix.append(ext2)
                self.__server_dict[str(ext2)] = []

            self.__server_dict[str(ext2)].append(str(item2)[:-3])
            prefix.sort()
            suffix.sort()
        #print('server_dict: ', self.__server_dict)
        if self.db_update_flag is False:
            ext_choice = self.u_i(
                ui_value,
                table=None,
                prefix=prefix,
                suffix=suffix)  # User Input
            f_t, e_c = ext_choice
        else:
            e_c = self.ext_choice
            f_t = self.file_choice
        for key in self.__server_dict.keys():
            temp_dict[str(key)] = []
            if str(key) not in str(f_t):
                del self.__server_dict[key]
                continue
            for ext in e_c:
                ext = str(ext)[:-3]
                print('ext: ', ext)
                for value in self.__server_dict[key]:
                    if str(ext) in str(value):
                        # print(value)
                        temp_dict[str(key)].append(value)
        self.__server_dict = {}
        #print('temp_dict: ', temp_dict)
        # __server_dict is a dictionary for our download files.
        self.__server_dict = temp_dict
        #print('server_dict: ', self.__server_dict)
        # The keys are the chosen file types and the values are the file names
        # with the selected extensions.
        return f_t, e_c

#  //TODO-ROB MODULE (DB):  Upload ftp files into a db file
    # //TODO-ROB give the user an option to do multiprocessing or not
    def db_upload(self, log_file):
        """Upload ftp files into a db file."""
        import json
        if self.database_flag == 'genbank':
            db_name = str(log_file).split('_', 1)
            db_name.pop(0)
            db_name = db_name[0]
            db_name = db_name.split('.')
            db_name = db_name.pop(0)

            files_path = db_name.replace('-', '/')
            print(files_path)
            files_path = where.NCBI_DATA + '/' + files_path + '/'
            ser_loc = where.DB
            t_count = 0
            home = os.getcwd()
            os.chdir(files_path)
            temp_list = []
            process = {}
            json_dict = {}
            with open(where.LOG + '/' + log_file, 'a+') as log_w:
                # //TODO-ROB Make better upload message for log file
                t = time.strftime("%m/%d/%Y")
                log_w.write(
                    '\n\n<%s>Below is a log of the databases created using a multiprocessing script:\n' %
                    t)
            print('Upload has started')
            for key, value in self.__server_dict.items():
                # //TODO-ROB Split __server_dict into multiple lists and send those list to multiple jobs

                print('key', key)
                print('value', value)

                process_count = 8
                div = len(self.__server_dict[key]) // (process_count - 1)
                nest_list = []
                count = 0
                db_count = 0
                for files in self.__server_dict[key]:
                    temp_list.append(files)
                    if len(
                            temp_list) % div == 0:  # 15 items per process for 16 processes
                        count += div
                        db_count += 1
                        nest_list.append(temp_list)
                        temp_list = []
                        if os.path.isfile(
                                ser_loc + ('/%s%s.%s.db' % (db_name, (db_count - 1), key))) is True:
                            print('DB file exists...  Continuing...')
                            continue
                        print('Copying Template BioSQL Database...  '
                              'This may take a few minutes...')
                        shutil.copy2(where.Templates + '/Template_BioSQL_DB.db',
                                     ser_loc + ('/%s%s.%s.db' % (db_name, (db_count - 1), key)))
                if (len(self.__server_dict[key]) - count) != 0:
                    db_count += 1
                    if os.path.isfile(ser_loc + ('/%s%s.%s.db' %
                                                 (db_name, (db_count - 1), key))) is True:
                        print('DB file exists...  Continuing...')
                    else:
                        shutil.copy2(where.Templates + '/Template_BioSQL_DB.db',
                                     ser_loc + ('/%s%s.%s.db' % (db_name, (db_count - 1), key)))
                    print('temp_list: ', temp_list)
                    nest_list.append(temp_list)
            os.chdir(where.DB)
            with open('temp_file.txt', mode='w', encoding='utf-8') as temp_file:
                json_dict["download_list"] = nest_list
                json_dict["db_name"] = db_name
                json_dict["key"] = key
                json_dict["log_file"] = str(log_file)
                json.dump(json_dict, temp_file)
        # with open(where.LOG + '/UPLOAD.log', 'a+') as lg_file:
        os.system('qsub UPLOAD.sh')

        print('Done submitting jobs')
        upload_flag = True
        while upload_flag == True:
            try:
                subprocess.check_output(['pidof', 'robupload'])
                upload_flag = False
            except subprocess.CalledProcessError:
                upload_flag = True
                time.sleep(30)
                print('Waiting for the upload to finish....')
        print('Uploading is complete.')
        os.remove('temp_file.txt')
        os.chdir(home)
        #//TODO-ROB Move the folder, and then move the db files into the directory


#  //TODO-ROB MODULE (DB): Fetch target GBK files from local NCBI datase.  Create a seperate class for this
#  //TODO-ROB MODULE (DB): Upload target GBK files into a new local database
#  //TODO-ROB MODULE (DB):  Update custom databases from the most recent donwloads
#------------------------------------------------------------------------------
    def config(self, config_file):
        """Update custom databases from the most recent downloads."""
        update_dict = {}
        self.__path_list = where.path_list_make(self.path, where.NCBI_DATA)
        os.sys.path.append(where.LOG)
        # Make a local dictionary from the log file called update_dict
        # and give it's values to the class variable __update_dict
        with open(where.LOG + '//' + config_file, 'r') as config_file:
            for line in log_file:
                if line.startswith('> update_dict'):
                    line = line.replace('> ', '')
                    print(line)
                    exec('%s' % line)
        self.__update_dict = update_dict

        # Configure all of the class variables with the dictionary
        for key, value in self.__update_dict.items():
            if isinstance(value, str):
                print(('ftp.' + str(key) + ' = ' + str(value)))
                # Initializes class string variables from the last archive
                exec(('self.%s = ' % key) + ('\'%s\'' % value))
                continue
            print('ftp.' + str(key) + ' = ' + str(value))
            # Initializes all other class variables from the last archive
            exec(('self.%s = ' % key) + ('%s' % value))

#------------------------------------------------------------------------------
# Update ****ALL**** NCBI files by using the most recent log files
    def update(self):
        if self.ftp_update_flag is True:
            for item in os.listdir(where.LOG):
                if str(self.__class_name).lower() in str(item).lower():
                    #//TODO-ROB Add JSON
                    # Class variables are configured from a log file
                    self.config(item)
                else:
                    continue
                self.ftp_update_flag = True
                os.system(
                    'rsync -vah --include \'*%s\' --exclude \'*\' ' %
                    ((self.__NCBI_RSYNC + '/' + self.path)))
                # Files are archived and the previous download is repeated from
                # log file
                self.ftp_download(self.path)
                input('ok')
        if self.db_update_flag is True:
            for log_file in self.__log_files:
                print('__log_files', log_file)
                self.config(log_file)
                self.ftp_extension(self.path)
                self.db_upload(log_file)
