import os
import sys
import zipfile
import logging as log
from pathlib import Path
from datetime import datetime as d
from Manager.utils.mana import ProjMana as PM
#import configparser
#from slacker import Slacker
#import argparse
#import textwrap
from Manager.utils.mana import UserMana


class DataMana(PM):

    def __init__(self, repo, user, project, research, research_type, home=os.getcwd(), data=None, new_data=False, new_project=False, **kwargs):
        super().__init__(repo=repo, user=user, project=project, research=research, research_type=research_type, home=home, new_project=new_project, **kwargs)
        self.data = data
        PM.


    def blast(self):
        print('blast folders')

    def genbank(self):
        print('genbank folders')

    def fasta(self):
        print('fasta foldrs')



# class LogIt(DataMana):
#     """LogIt makes logging easier by creating easy loggers."""
#     """Class for improving the ease of use of the logging module."""
#     def __init__(self, logfile=None, logname=None):
#         """Initialize the logger format based on system platform."""
#         # Set the different formats
#         """Set the different formats based on user's platform."""
#         if sys.platform == 'win32':
#             self.archive_format = '%m-%d-%Y_%I-%M-%p'
#             pass
#         elif sys.platform == 'linux':
#             self.archive_format = '%m-%d-%Y@%I:%M:%S-%p'
#             pass
#
#         self.date_format = '%a %b %d at %I:%M:%S %p %Y'  # Used to add as a date
#         self.log_format = '%(name)s - [%(levelname)-2s]: %(message)s'
#         # self.slack = self.slack_config()
#         self.basic = self.generic_logger(logfile, logname, log.DEBUG, self.log_format)
#
#     def _get_file(self, filename):
#         """Create a log file."""
#         base, extension = filename.split('.')
#         file = base + str(d.now().strftime(self.archive_format)) + extension
#         path = os.getcwd() + file
#         return path
#
#     def generic_logger(self, filename, logname, level, fmt, slack=False):
#         """Create a generic logger."""
#         file_path = self._get_file(filename)
#         log.basicConfig(level=level,
#                         format=fmt,
#                         filename=file_path)
#         generic_logger = log.getLogger(logname)
#         if slack is False:
#             return generic_logger
#             #        else:
#             #            slack_logger = log.getLogger('SLACK')
#
#             # ******************************************SLACK****************************************** #
#             # ******************************************SLACK****************************************** #
#             # ******************************************SLACK****************************************** #
#
#             # @staticmethod
#             # def slack_config():
#             #     config = configparser.ConfigParser()
#             #     config.read('bin/orthologs.ini')
#             #     apikey = config['APIKEYS']['slack']
#             #     slack = Slacker(apikey)
#             #     return slack
#             #
#             # # Definition for uploading images
#             # def upload_img(self, channel, imgfile):
#             #     self.slack.files.upload(imgfile, channel=channel)
#             #
#             # # Definition for uploading files
#             # def upload_file(self, channel, file):
#             #     self.slack.files.upload(file, channel=channel)
#             #
#             # # Definition for posting messages
#             # def message_slack(self, channel, message, username):
#             #     self.slack.chat.post_message(channel, message, username, as_user=True)
#
#             # ******************************************SLACK****************************************** #
#             # ******************************************SLACK****************************************** #
#             # ******************************************SLACK****************************************** #
#
#
# class ZipUtils(DataMana):
#     """The ZipUtil class allows easy compression/zipping of file folders.
#     Inspired by http://stackoverflow.com/a/670635/7351746
#     """
#
#     def __init__(self, zip_filename, zip_path):
#         """Initialize the input files and path.
#
#         :param comp_filename (string):  This is the name of the compressed file that will be generated (eg 'test.zip')
#         :param zip_path: This is the absolute path of the directory (or file) to be zipped.
#         :returns:  A zip file that is created inside of the zip_path.  The path string is returned.
#         """
#         self.zip_filename = zip_filename
#         self.zip_path = zip_path
#         self.ignore_parts = Path(zip_path).parent.parts
#
#     def to_zip(self):
#         """Zip a folder."""
#         comp_path = os.path.join(self.zip_path, self.zip_filename)
#         zip_handle = zipfile.ZipFile(comp_path, 'w', zipfile.ZIP_DEFLATED)
#         if os.path.isfile(self.zip_path):
#             zip_handle.write(self.zip_path)
#         else:
#             print('skipped')
#             self.add_folder_to_zip(zip_handle, self.zip_path)
#         zip_handle.close()
#         return comp_path
#
#     def add_folder_to_zip(self, zip_handle, folder):
#         """Not meant to be used explicitly.  Use to_zip.
#
#         :param zip_handle: An initialized zipfile.ZipFile handle.
#         :param folder: A path that represents an entire folder to be zipped.
#         :return: Recursively zips nested directories.
#         """
#         for file in os.listdir(folder):
#             full_path = os.path.join(folder, file)
#             rel_path = Path(full_path)
#             rel_path = rel_path.relative_to(Path(self.zip_path))
#             if os.path.isfile(full_path):
#                 if str(file) == str(self.comp_filename):
#                     continue
#                 print('File added: ' + str(full_path))
#                 zip_handle.write(full_path, rel_path)
#             elif os.path.isdir(full_path):
#                 if str(file) in self.ignore_parts:
#                     continue
#                 print('Entering folder: ' + str(full_path))
#                 self.add_folder_to_zip(zip_handle, full_path)
