"""Main logging class to make logging easier."""

import logging as log
import os
from datetime import datetime as d
import sys
#import configparser
#from slacker import Slacker
#import argparse
#import textwrap

class LogIt(object):
    """LogIt makes logging easier by creating easy loggers."""
    def __init__(self, logfile=None, logname=None):
        """Initialize the logger format based on system platform."""
        # Set the different formats
        if sys.platform == 'win32':
            self.archive_format = '%m-%d-%Y_%I-%M-%p'
            pass
        elif sys.platform == 'linux':
            self.archive_format = '%m-%d-%Y@%I:%M:%S-%p'
            pass

        self.date_format = '%a %b %d at %I:%M:%S %p %Y'  # Used to add as a date
        self.log_format = '%(name)s - [%(levelname)-2s]: %(message)s'
        # self.slack = self.slack_config()
        self.basic = self.generic_logger(
            logfile, logname, log.DEBUG, self.log_format)

    def _get_file(self, filename):
        """Create a log file."""
        base, extension = filename.split('.')
        file = base + str(d.now().strftime(self.archive_format)) + extension
        path = os.getcwd() + file
        return path

    def generic_logger(self, filename, logname, level, fmt, slack=False):
        """Create a generic logger."""
        file_path = self._get_file(filename)
        log.basicConfig(level=level,
                        format=fmt,
                        filename=file_path)
        generic_logger = log.getLogger(logname)
        if slack is False:
            return generic_logger
#        else:
#            slack_logger = log.getLogger('SLACK')

    # @staticmethod
    # def slack_config():
    #     config = configparser.ConfigParser()
    #     config.read('bin/orthologs.ini')
    #     apikey = config['APIKEYS']['slack']
    #     slack = Slacker(apikey)
    #     return slack
    #
    # # Definition for uploading images
    # def upload_img(self, channel, imgfile):
    #     self.slack.files.upload(imgfile, channel=channel)
    #
    # # Definition for uploading files
    # def upload_file(self, channel, file):
    #     self.slack.files.upload(file, channel=channel)
    #
    # # Definition for posting messages
    # def message_slack(self, channel, message, username):
    #     self.slack.chat.post_message(channel, message, username, as_user=True)