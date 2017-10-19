"""Collection of tools for using PBS, a job scheduler for high-performance
computing environments. The command is usually `qsub <options>` on most
systems.
"""
import random
import string
from string import Template
from subprocess import run, CalledProcessError, PIPE
from datetime import datetime as d
import platform
import getpass
import os
import contextlib


class Qsubutils:
    """Create a pbs job and submit it using qsub.

    This class also provides functionality for creating multiple pbs jobs that
    by creating chunks of lists for each python script and job.
    """
    def __init__(self, default=True):
        """UNLESS A WINDOWS MACHINE HAS PBS. IT LIKELY WONT"""
        if "windows" in platform.system().lower():
            raise ImportError("QsubUtils is only supported on linux/osx.")

        self.default = default

    def basejobids(self):
        """"Create base job attributes."""
        base_id = self.randomid()
        self.base_id = base_id
        base = "submit_{0}".format(self.base_id)

        return base_id, base

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
    def randomid(length=5):
        """Generate a random ID of 5 characters to append to qsub job name."""
        return ''.join(random.sample(
            string.ascii_letters + string.digits, length))

    @staticmethod
    def writecodefile(filename, code, language):
        """Create a python file and write the code to it."""
        if language == 'python':
            with open(filename + '.py', 'w') as pyfile:
                pyfile.write(code)
                pyfile.close()

        elif language == 'bash':
            with open(filename + '.sh', 'w') as bashfile:
                bashfile.write(code)
                bashfile.close()

        elif language == 'R' or 'r':
            with open(filename + '.R', 'w') as rfile:
                rfile.write(code)
                rfile.close()
        else:
            raise NotImplementedError('%s is unsupported.' % language)

    def _checkjobstatus(self):
        raise NotImplementedError
#        with contextlib.suppress(OSError):
#            cmd = 'qstat'  # this is the command
#            cmd_status = run([cmd], shell=True)  # must = TRUE
#            if cmd_status == 0:  # Command was successful.
#                print('Job submitted.')
#            else:  # Unsuccessful. Stdout will be '1'.
#                print("PBS job not submitted.")

    @classmethod
    def _cleanup(cls, base):
        os.remove(base + '.pbs')
        os.remove(base + '.py')

    def pbs_dict(self, default=True, author='', email='', description='', project='', select='1',
                 jobname='', script_name='', logname='', outfile='', errfile='', python='python3',
                 gigabytes='6gb', cput='75:00:00', walltime='75:00:00'):
        """Add PBS script attributes."""

        # TODO simplify this function. Add configuration script.

        # Date format. Used to add as a date
        format1 = '%a %b %d %I:%M:%S %p %Y'

        if self.default == default:
            author = getpass.getuser().upper()
            email = 'n/a'
            description = 'This is a default pbs job.'
            project = 'Datasnakes-Orthologs'
            select = '3'
            jobname = 'orthoevol'
            script_name = script_name
            outfile = outfile
            errfile = errfile
            logname = 'ortho-evol.{}'.format(self.randomid(length=3))

        # Fill in the pbs script attributes
        job_attributes = {
            'author': author,
            'description': description,
            'date': d.now().strftime(format1),
            'PBS_JOBID': '${PBS_JOBID}',
            'PBS_O_WORKDIR': '${PBS_O_WORKDIR}',
            'proj_name': project,
            'select': int(select),
            'memgb': gigabytes,
            'cput': cput,
            'wt': walltime,
            'job_name': jobname,
            'outfile': outfile,
            'errfile': errfile,
            'script': script_name,
            'log_name': logname.format(self.randomid()),
            'cmd': python + " " + script_name + ".py",
            'email': email
            }

        return job_attributes

    def submitjob(self, code, language='python', default=True, prefix=None):
        """Creates and submit a qsub job. Also uses python code."""
        # TIP If python is in your environment as only 'python' update that.
        # TODO-SDH add a default parameter or custom parameter
        # If default, a python file will be created from code that is used.
        if self.default == default and language == 'python':
            baseid, base = self.basejobids()
            if prefix is not None:
                base = prefix + '_' + base
            self.writecodefile(filename=base, code=code, language=language)
            outfile = 'orthoevol_{}.out'.format(baseid)
            errfile = 'orthoevol_{}.err'.format(baseid)
            # Create the pbs script from the template or dict
            pbstemp = self.import_temp('temp.pbs')

            script_name = base.format(baseid)

            attributes = self.pbs_dict(outfile=outfile,
                                       errfile=errfile,
                                       script_name=script_name)

            with open(base + '.pbs', 'w') as pbsfile:
                pbsfile.write(pbstemp.substitute(attributes))
                pbsfile.close()
        else:
            raise NotImplementedError('Custom qsub jobs are forbidden.')
            # TODO Improve custom job creation
#            pbstemp = self.import_temp('temp.pbs')
#            with open(base + '.pbs', 'w') as pbsfile:
#                pbsfile.write(pbstemp.substitute(self.pbs_dict))
#                pbsfile.close()

        with contextlib.suppress(CalledProcessError):
            cmd = ['qsub ' + base + '.pbs']  # this is the command
            # Shell MUST be True
            cmd_status = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

            if cmd_status.returncode == 0:  # Command was successful.
                print('Job submitted.\n')
                # TODO add a check to for job errors or check for error file.

            else:  # Unsuccessful. Stdout will be '1'
                print("PBS job not submitted.")
                self._cleanup(base=base)