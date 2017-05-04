# -*- coding: utf-8 -*-
# Modules Used
import os
import random
import string
#import tempfile
import subprocess
from string import Template
from datetime import datetime as d

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
    file_temp = open(filename, 'r')
    file_str = file_temp.read()
    return file_str

#------------------------------------------------------------------------------
def RandomId(length=5):
    return ''.join(random.sample(string.ascii_letters + string.digits, length))

#------------------------------------------------------------------------------
def SubmitPythonCode(code, pbstemp, author, jobname="job", cleanup=True, prefix="", slots=1):
    """ This function creates and submits qsub jobs."""

    format1 = '%a %b %d %I:%M:%S %p %Y'  # Used to add as a date

    # Create a base name for the jobs and python scripts
    base = prefix + "submit_{0}".format(RandomId())
    open(base + '.py', 'w').write(code)
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

    open(base + '.pbs', 'w').write(pbstemp.substitute(pbs_dict))
    try:
        cmd = 'qsub  ' + base + '.pbs'
        cmd_status = subprocess.call([cmd], shell=True)
        if cmd_status == 0:  # Command was successful.
            pass  # Continue
        else:  # Unsuccessful. Stdout will be '1'.
            print("PBS job not submitted.")
    finally:
        if cleanup:
            os.remove(base + '.qsub')
            os.remove(base + '.py')

#------------------------------------------------------------------------------

