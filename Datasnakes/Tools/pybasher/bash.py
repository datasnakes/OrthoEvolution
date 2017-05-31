"""A user-friendly Bash module for Python."""
# TODO-SDH It may be helpful to use psutil for pybasher.
# TODO-SDH Finish pybasher ASAP
# TODO-SDH Look at some examples for pybasher.
# TODO-SDH Update the README for pybasher
import platform
import sys


SUBPROCESS_HAS_TIMEOUT = True
if "windows" in platform.system().lower():
    raise ImportError(
            "sh %s is currently only supported on linux and osx.")
elif sys.version_info < (3, 0):
    try:
        from subprocess32 import PIPE, Popen
    except ImportError:
        # You haven't got subprocess32 installed. If you're running 2.X this
        # will mean you don't have access to things like timeout
        SUBPROCESS_HAS_TIMEOUT = False

from subprocess import PIPE, Popen, call
#import os
#import configparser
# TODO-SDH use a config file to load/use a list or group of common commands.


class PyBasher(object):
    # !!! Only for linux
    def __init__(self, *args, **kwargs):
        """Initialize the call as well as standard error and output."""
        # TODO-SDH Test if this is working.
        self.process = None
        self.stdout = None
        self.pybasher(*args, **kwargs)
        
    def runcmd(self, cmd, env=None, stdout=PIPE, stderr=PIPE, timeout=None, sync=True):
        self.process = Popen(cmd, shell=True, stdout=stdout, stdin=PIPE, stderr=stderr, env=env)
        if sync:
            self.sync(timeout)
        return self