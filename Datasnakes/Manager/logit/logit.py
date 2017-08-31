"""Main logging class to make logging easier."""
import logzero
from logzero import logger as log
import os
from datetime import datetime as d
import sys
import configparser
#from slacker import Slacker



class LogIt(object):
    """LogIt makes logging easier by creating easy loggers."""
    def __init__(self, logfile=None, logname=None):
        """Initialize the logger format based on system platform."""
        # Set the different formats based on user's platform
        if sys.platform == 'win32':
            self.archive_format = '%m-%d-%Y_%I-%M-%p'
        elif sys.platform == 'linux':
            self.archive_format = '%m-%d-%Y@%I:%M:%S-%p'

        self.date_format = '%a %b %d at %I:%M:%S %p %Y'  # Used to add date
        self.log_format = '%(name)s - [%(levelname)-2s]: %(message)s'
        # self.slack = self.slack_config()
        self.basic = self.generic_logger(
            logfile, logname, log.DEBUG, self.log_format)

