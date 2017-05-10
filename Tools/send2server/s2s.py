# -*- coding: utf-8 -*-
"""
File Name:
Description: s2s sets up sending files to servers via public SSH keys.

Author: S. Hutchins
Date Created: Tue Apr 11 12:31:17 2017
Project Name: Orthologs Project
"""

#import os
#import logging as log
#import pandas as pd
#from datetime import datetime as d
##import zipfile
#import pexpect
import subprocess

#------------------------------------------------------------------------------
class S2S(object):
    """S2S (Send 2 Server) is designed for use with a public ssh key."""
    def __init__(username, server_address):
        address = server_address
        user = username
        sendto = user + "@" + address + ":"
        return sendto

    def scpto(self, file, destpath):
        cmd = "scp " + file + " " + self.sendto + destpath
        status = subprocess.call([cmd], shell=True)
        if status == 0:  # Command was successful.
            print("%s file sent." % file)
            pass  # Continue
        else:  # Unsuccessful. Stdout will be '1'.
            print("%s file not sent." % file)
