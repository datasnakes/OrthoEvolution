# -*- coding: utf-8 -*-
"""
File Name: SubmitQsubJob.py
Description:

Author: S. Hutchins
Date Created: Thu May  4 16:18:08 2017
Project Name: Orthologs Project
"""

# Modules Used
import os
import random
import string
#import tempfile  # This MAY be added later.
import subprocess
from string import Template
from datetime import datetime as d
import sys

#------------------------------------------------------------------------------
class CreateJob(object):
    def __init__():
        """UNLESS A WINDOWS MACHINE HAS PBS (IT SHOULDNT)"""
        if sys.platform == 'win32' or 'win64':
            sys.exit("This module is strictly for use on Linux at the moment.")
#------------------------------------------------------------------------------
    def ImportTemp(filepath):
        """Import the script or file that you need a template of and that has temp
        strings."""
        file_temp = open(filepath, 'r')
        file_str = file_temp.read()
        file_temp.close()

        file_temp = Template(file_str)
        return file_temp

#------------------------------------------------------------------------------
    def File2Str(filename):
        """ Turn the contents of a file (python file) into a string."""
        file_temp = open(filename, 'r')
        file_str = file_temp.read()
        return file_str

#------------------------------------------------------------------------------
    def RandomId(self, length=5):
        """Generate a random ID of 5 characters to append to qsub job name."""
        return ''.join(random.sample(string.ascii_letters + string.digits, length))

#------------------------------------------------------------------------------
    def SubmitPythonCode(self, code, pbstemp, author, jobname="job", cleanup=True, prefix="", slots=1):
        """ This function creates and submits qsub jobs."""

        format1 = '%a %b %d %I:%M:%S %p %Y'  # Used to add as a date

        # Create a base name for the jobs and python scripts
        base = prefix + "submit_{0}".format(self.RandomId())

        # Create a python file and write the code to it
        open(base + '.py', 'w').write(code)

        # This is the "command" for running the python script. If python is in
        # your environment as only 'python' update that here.
        script = "python3 " + base + ".py"

        pbs_dict = {
                'author': author,
                'description': 'do things',
                'date': d.now().strftime(format1),
                'PBS_JOBID': '${PBS_JOBID}',
                'PBS_O_WORKDIR': '${PBS_O_WORKDIR}',
                'proj_name': 'Orthologs Project',
                'select': '4',
                'memgb': '8gb',
                'cput': '30:00:00',
                'wt': '30:00:00',
                'job_name': jobname,
                'script': script,
                'log_name': base,
                'cmd': ''
                }
        # Create the pbs script from the template/dict
        open(base + '.pbs', 'w').write(pbstemp.substitute(pbs_dict))

        # Submit the qsub job using subprocess
        try:
            cmd = 'qsub  ' + base + '.pbs'  # this is the command
            cmd_status = subprocess.call([cmd], shell=True)  #  Shell must be TRUE
            if cmd_status == 0:  # Command was successful.
                pass  # Continue
            else:  # Unsuccessful. Stdout will be '1'.
                print("PBS job not submitted.")
        finally:
            if cleanup:  # When finished, remove the qsub files and python files.
                os.remove(base + '.qsub')
                os.remove(base + '.py')

#------------------------------------------------------------------------------

