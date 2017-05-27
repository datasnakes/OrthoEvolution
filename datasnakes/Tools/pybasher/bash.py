"""A user-friendly Bash module for Python."""
# TODO-SDH It may be helpful to use psutil for pybasher.
# TODO-SDH Finish pybasher ASAP
# TODO-SDH Look at some examples for pybasher.
# TODO-SDH Update the README for pybasher
from subprocess import PIPE, Popen
import sys
#import os
#import configparser
# TODO-SDH use a config file to load/use a list or group of common commands.


class PyBasher(object):
    # !!! Only for linux
    def __init__(self, cmd):
        """Initialize the call as well as standard error and output."""
        # TODO-SDH Test if this is working.
        c = Popen([cmd], shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = c.communicate()
        if stderr:
            sys.stderr
        else:
            sys.stdout          
#        return c
#
#        __COMMANDS__ = [
#                ]
#        commands = __COMMANDS__
#        return commands
#    
#    def listcommands(configfile='bash.cfg'):
#        """Use a bash configuration file of common commands."""
#        # IDEA not sure if this is the best option
