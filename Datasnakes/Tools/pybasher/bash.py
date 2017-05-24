# -*- coding: utf-8 -*-
"""
File Name:
Description:

Author: S. Hutchins
Date Created: Mon May  8 11:09:04 2017
Project Name: Orthologs Project
"""
from subprocess import call, check_output
import os

#------------------------------------------------------------------------------


class CleanUp(object):
    def __init__(self, cmd):
        c = call([cmd], shell=True)
        return c

    def rmdir(self, cmd="rm -r " + pathname, pathname):
        self.c
        return print("%s was removed." % pathname)
