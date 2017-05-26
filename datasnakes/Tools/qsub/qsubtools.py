"""Collection of tools for using PBS, a job scheduler for high-performance
computing environments. The command is usually `qsub <options>` on most
systems.
"""
import os
import random
import string
import subprocess
from string import Template
from datetime import datetime as d
import sys


class CreateJob(object):
    """Create a pbs job and submit it using qsub.

    This class also provides functionality for creating multiple pbs jobs that
    by creating chunks of lists for each python script and job.
    """
    def __init__():
        """UNLESS A WINDOWS MACHINE HAS PBS IT LIKELY WONT"""
        # XXX This class is optimized for Linux!
        if sys.platform == 'win32' or 'win64':
            sys.exit("This module is strictly for use on Linux at the moment.")

    def import_temp(filepath):
        """Import the script or file that you need a template of and that has
        temp strings.
        """
        file_temp = open(filepath, 'r')
        file_str = file_temp.read()
        file_temp.close()

        file_temp = Template(file_str)
        return file_temp

    def file2str(filename):
        """Turn the contents of a file (python file) into a string."""
        file_temp = open(filename, 'r')
        file_str = file_temp.read()
        return file_str

    def randomid(self, length=5):
        """Generate a random ID of 5 characters to append to qsub job name."""
        return ''.join(random.sample(
            string.ascii_letters + string.digits, length))

    def pbs_dict(self, author, email, description, project, select, jobname,
                 script_name, logname, python='python3', gb='6gb',
                 cput='75:00:00', walltime='75:00:00'):
        """Add PBS script attributes."""
        # Date format
        format1 = '%a %b %d %I:%M:%S %p %Y'  # Used to add as a date
        # Fillin the pbs script attributes
        job_attributes = {
            'author': author,
            'description': description,
            'date': d.now().strftime(format1),
            'PBS_JOBID': '${PBS_JOBID}',
            'PBS_O_WORKDIR': '${PBS_O_WORKDIR}',
            'proj_name': project,
            'select': int(select),
            'memgb': gb,
            'cput': cput,
            'wt': walltime,
            'job_name': jobname,
            'script': script_name,
            'log_name': logname,
            'cmd': python + " " + script_name + ".py",
            'email': email
            }

        return job_attributes

    def submitpythoncode(self, code, cleanup=True, prefix="", slots=1):
        """Creates and submit a qsub job. Also use python code."""

        # Create a base name for the jobs and python scripts
        base = prefix + "submit_{0}".format(self.randomid())

        # Create a python file and write the code to it
        with open(base + '.py', 'w') as pyfile:
            pyfile.write(code)
            pyfile.close()

        # This is the "command" for running the python script.
        # TIP If python is in your environment as only 'python' update that.
        # TODO-SDH add a default parameter or custom parameter
        # Custom will allow user to use whatever python script they want
        # If default, python file will be created from code that is used.
        script = python + " " + base + ".py"

        # Create the pbs script from the template or dict
        with open(base + '.pbs', 'w') as pbsfile:
            pbsfile.write(pbstemp.substitute(self.pbs_dict))
            pbsfile.close()

        # Submit the qsub job using subprocess
        try:
            cmd = 'qsub  ' + base + '.pbs'  # this is the command
            print(cmd)
#            cmd_status = subprocess.call([cmd], shell=True)  # must = TRUE
#            if cmd_status == 0:  # Command was successful.
#                pass  # Continue
#            else:  # Unsuccessful. Stdout will be '1'.
#                print("PBS job not submitted.")
        finally:
            print('finally')
#            if cleanup:  # When finished, remove the qsub files & python files.
#                os.remove(base + '.qsub')
#                os.remove(base + '.py')

    def multiplejobs(self):
        """Create multiple jobs & scripts for each job to run based on
        splitting a list into chunks.
        """
        print(self.__doc__)
#        # Modules used
#        from QsubTools import SubmitPythonCode, ImportTemp
#        from multiprocess_functions import genes2align
#
#        #------------------------------------------------------------------------------
#        # Import the pbs and script templates
#        pbs_temp = ImportTemp(filepath='templates/temp.pbs')
#
#        script_temp = """from multiprocess_functions import main, clustal
#        main(geneslist={}, function=clustal)
#            """
#
#        list_chunks = genes2align()
#
#        for k, v in list_chunks.items():
#            # Format the script template with the list of genes to align or run PAML on
#            # Submit the python code and create qsub jobs
#            SubmitPythonCode(code=script_temp.format(v), pbstemp=pbs_temp, author="SDH", jobname="karg-clust")
#
#        # The next major part to add is creation of a bash script or waiting for
#        # the genes to align in order to start paml
