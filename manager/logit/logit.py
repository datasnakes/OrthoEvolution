# -*- coding: utf-8 -*-
"""
Date created: Sat Apr 22 20:45:17 2017
Author: S. Hutchins

Description: Main logging class

"""
import logging as log
import os
from datetime import datetime as d
import sys


class LogIt(object):

    def __init__(self, logfile=None, logname=None):
        #
        # # Variables for logging formats
        # message_format = '%(name)s - [%(levelname)-2s]: %(message)s'
        #
        #
        if sys.platform == 'win32':
            self.archive_format = '%m-%d-%Y_%I-%M-%p'
            pass
        elif sys.platform == 'linux':
            self.archive_format = '%m-%d-%Y@%I:%M:%S-%p'
            pass
        #
        # # Would be a great idea to add a log directory setup if that directory
        # # doesn't exist.
        #
        # log.basicConfig(level=log.DEBUG,
        #                 format=message_format,
        #                 filename="logs/" + logfile + "_%s.log" % str(d.now().strftime(logappend_format)))
        # self.log = log.getLogger(logname)
        # Logging variables
        self.date_format = '%a %b %d at %I:%M:%S %p %Y'  # Used to add as a date
        # self.archive_format = '%m-%d-%Y@%I:%M:%S-%p'  # Used to append to archives
        self.log_format = '%(name)s - [%(levelname)-2s]: %(message)s'

    def blastn(self):
        log.basicConfig(level=log.DEBUG,
                        format=self.log_format,
                        filename="%s/BLAST_%s.log" % (os.getcwd(), str(d.now().strftime(self.archive_format))))
        self.blast_log = log.getLogger('Blastn')
        return self.blast_log

    def post_blast(self):
        log.basicConfig(level=log.INFO,
                        format=self.log_format,
                        filename="%s/Post_BLAST_Analysis_%s.log" % (os.getcwd(), str(d.now().strftime(self.archive_format))))
        self.post_blast_log = log.getLogger('PostBlast')
        return self.post_blast_log

    def scriptinfo(self):
        # Write basic information to the log
        date_format = '%a %b %d at %I:%M:%S %p %Y'  # Used to add as a date
        # Write basic information to the log
        self.log.info("------------------------------------------------------------------")
        self.log.info("The script name is %s" % os.path.basename(sys.argv[0]))
        self.log.info("The script began on %s" % str(d.now().strftime(date_format)))
        self.log.info("------------------------------------------------------------------")



