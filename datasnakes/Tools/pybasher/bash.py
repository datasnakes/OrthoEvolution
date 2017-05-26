"""A user-friendly Bash module for Python."""
# TODO-SDH It may be helpful to use psutil for pybasher.
# TODO-SDH Finish pybasher ASAP
# TODO-SDH Look at some examples for pybasher.
# TODO-SDH Update the README for pybasher
from subprocess import call, check_output
import os


class PyBasher(object):
    def __init__(self, cmd):
        # TODO-SDH Test if this is working.
        c = call([cmd], shell=True)
        return c

    def rmdir(self, cmd="rm -r " + pathname, pathname):
        self.c
        return print("%s was removed." % pathname)
