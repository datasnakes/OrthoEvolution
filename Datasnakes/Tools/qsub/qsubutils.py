"""Collection of tools for using PBS, a job scheduler for high-performance
computing environments. The command is usually `qsub <options>` on most
systems.
"""
import sys
import random
import string
from string import Template
import subprocess
from datetime import datetime as d
import platform
import getpass
import os
from Datasnakes.Tools.multiprocess import genes2analyze


class QsubUtils:
    """Create a pbs job and submit it using qsub.

    This class also provides functionality for creating multiple pbs jobs that
    by creating chunks of lists for each python script and job.
    """
    def __init__(self):
        """UNLESS A WINDOWS MACHINE HAS PBS. IT LIKELY WONT"""
        if "windows" in platform.system().lower():
            raise ImportError("QsubUtils is only supported on linux/osx.")

        default = True  # Default job is TRUE
        self.default = default

    @classmethod
    def import_temp(cls, filepath):
        """Import the script or file that you need a template of and that has
        temp strings.
        """
        file_temp = open(filepath, 'r')
        file_str = file_temp.read()
        file_temp.close()

        file_temp = Template(file_str)
        return file_temp

    @classmethod
    def file2str(cls, filename):
        """Turn the contents of a file (python file) into a string."""
        file_temp = open(filename, 'r')
        file_str = file_temp.read()
        return file_str

    @staticmethod
    def _randomid(length=5):
        """Generate a random ID of 5 characters to append to qsub job name."""
        return ''.join(random.sample(
            string.ascii_letters + string.digits, length))

    def pbs_dict(self, default, author='', email='', description='', project='', select='1',
                 jobname='', script_name='', logname='', python='python3',
                 gb='6gb', cput='75:00:00', walltime='75:00:00'):
        """Add PBS script attributes."""
        format1 = '%a %b %d %I:%M:%S %p %Y'  # Date format. Used to add as a date

        if self.default == default:
            author = getpass.getuser().upper()
            email = 'n/a'
            description = 'This is a default pbs job.'
            project = 'Datasnakes-Orthologs'
            select = '3'
            jobname = 'ds-ortho'
            script_name = 'datasnakes'
            logname = 'datasnakes'.format(self._randomid(length=3))

        # Fill in the pbs script attributes
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
            'log_name': logname.format(self._randomid()),
            'cmd': python + " " + script_name + ".py",
            'email': email
            }

        return job_attributes

    def submitpythoncode(self, code, default=True, cleanup=True, prefix=""):
        """Creates and submit a qsub job. Also uses python code."""

        # Create a base name for the jobs and python scripts
        base = prefix + "submit_{0}".format(self._randomid())

        # Create a python file and write the code to it
        with open(base + '.py', 'w') as pyfile:
            pyfile.write(code)
            pyfile.close()
        # TIP If python is in your environment as only 'python' update that.
        # TODO-SDH add a default parameter or custom parameter
        # If default, python file will be created from code that is used.
        if self.default == default:
            # Create the pbs script from the template or dict
            pbstemp = self.import_temp('temp.pbs')
            with open(base + '.pbs', 'w') as pbsfile:
                pbsfile.write(pbstemp.substitute(self.pbs_dict(self.default)))
                pbsfile.close()
        else:
            pbstemp = self.import_temp('temp.pbs')
            with open(base + '.pbs', 'w') as pbsfile:
                pbsfile.write(pbstemp.substitute(self.pbs_dict(default=False)))
                pbsfile.close()

        userprompt = input('Are you sure you want to submit this job? [Y/N]')
        if userprompt != 'Y' or 'y':
            # Delete the jobs that were just created if job not performed.
            os.remove(base + '.pbs')
            os.remove(base + '.py')
            sys.exit('You do not want to submit this job.')  # Exit.

        try:
            cmd = 'qsub  ' + base + '.pbs'  # this is the command
            cmd_status = subprocess.call([cmd], shell=True)  # must = TRUE
            if cmd_status == 0:  # Command was successful.
                pass  # Continue
            else:  # Unsuccessful. Stdout will be '1'.
                print("PBS job not submitted.")
        finally:
            print('hi')
            if cleanup:  # When finished, remove the qsub files & python files.
                os.remove(base + '.qsub')
                os.remove(base + '.py')


class MultiJobber(QsubUtils):
    """Create multiple jobs & scripts for each job to run based on
    splitting a list into chunks.
    """
    # TODO-SDH This needs testing.
    # TODO-SDH Create a more simplified process/algorithm.
    def __init__(self, default=True, cleanup=True, prefix=""):
        super().__init__(self)
        self.code = super().import_temp('temp.py')
        create = super().submitpythoncode(self, code=self.code,
                                          default=True,
                                          cleanup=True,
                                          prefix="")

        self.clustal_chunks, self.paml_chunks = genes2analyze()

        try:
            for key, value in self.clustal_chunks.items():
                create
        except:
            raise Exception('Error in clustal multijobber.')
            sys.exit('Will not proceed to paml.')
        try:
            for key, value in self.paml_chunks.items():
                create
        except:
            raise Exception('Error in pamlmultijobber.')
            sys.exit('Will not proceed')

        finally:
            print('Clustal and Paml multijobber complete.')
