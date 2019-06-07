import getpass
import string
import random
import os
import subprocess as sp
from pathlib import Path
from collections import OrderedDict
from pkg_resources import resource_filename
from datetime import datetime as d
from time import sleep
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Manager.config import templates
from OrthoEvol.utilities import FullUtilities
from OrthoEvol.Tools.sge import Qstat


class BaseQsub(object):
    """Base class for PBS jobs."""

    def __init__(self, pbs_script=None):

        self.pbs_log = LogIt().default(logname="PBS JOB", logfile=None)
        self.pbs_utils = FullUtilities()
        self.pbs_script = pbs_script

    def submit_pbs_script(self):
        """Submit a job using qsub.

        :param cleanup: (Default value = False)
        :param wait: (Default value = True)
        """
        try:
            cmd = ['qsub ' + str(self.pbs_script)]  # this is the command
            # Shell MUST be True
            cmd_status = self.pbs_utils.system_cmd(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, check=True)
        except sp.CalledProcessError as err:
            self.pbs_log.error(err.stderr.decode('utf-8'))
        else:
            if cmd_status.returncode == 0:  # Command was successful.
                # The cmd_status has stdout that must be decoded.
                # When a qsub job is submitted, the stdout is the job id.
                submitted_jobid = cmd_status.stdout.decode('utf-8')
                self.pbs_log.info(str(self.pbs_script) + ' was submitted.')
                self.pbs_log.info('Your job id is: %s' % submitted_jobid)

            else:  # Unsuccessful. Stdout will be '1'
                self.pbs_log.error('PBS job not submitted.')


class Qsub(BaseQsub):

    def __init__(self, author=getpass.getuser(), project_name="OrthoEvol", description="This is a basic pbs job",
                 date_format='%a %b %d %I:%M:%S %p %Y', chunk_resources=None, cput='72:00:00', walltime='48:00:00',
                 job_name=None, pbs_work_dir=None, email=None, directive_list=None, primary_pbs_cmd=None,
                 pbs_command_list=None, email_command=None, activate_script=None):

        self.base_name = job_name
        self.base_id, self.job_name = self.get_base_job_name()
        self.pbs_script = Path(self.pbs_work_dir) / self.job_name + '.pbs'
        self.python_script = Path(self.pbs_work_dir) / self.job_name + '.py'
        super().__init__(pbs_script=self.pbs_script)

        self.pbs_template = resource_filename(templates.__name__, "temp.pbs")

        # PBS - header info
        self.author = author
        self.project_name = project_name
        self.description = description
        self.current_date = d.now().strftime(date_format)
        # PBS - directives/attributes
        self.resource_str = self.get_resource_string(chunk_resources=chunk_resources)
        self.cputime = cput
        self.walltime = walltime
        self.directive_list = directive_list
        # PBS - commands
        self.pbs_work_dir = pbs_work_dir
        self.activate_script = activate_script  # for python virtual environment activation
        self.primary_pbs_cmd = primary_pbs_cmd
        self.pbs_command_list = pbs_command_list
        self.email = email
        if email_command is None:
            if email is not None:
                self.email_command = "mail -s \"%s.py script completed\" %s <<< 'Check your output'" % \
                                     (self.job_name, self.email)
            else:
                self.email_command = None
        else:
            self.email_command = email_command

    def get_resource_string(self, chunk_resources=None):
        if chunk_resources is None:
            chunk_resources = OrderedDict({
                "select": 1,
                "ncpus": 1,
                "memgb": "6gb"
            })

        resource_list = []
        for k, v in chunk_resources.items():
            if v is not None:
                resource_list.append("%s=%s" % (k, v))
        resource_str = ":".join(resource_list)
        return resource_str

    def get_base_job_name(self, length=5):
        base_id = ''.join(random.sample(string.ascii_letters + string.digits, length))
        job_name = self.base_name + "_%s" % base_id

        return base_id, job_name

    def _cleanup(self, jobname):
        """Clean up job scripts.

        :param jobname: The name of the job being run or to be run.
        """

        self.pbs_log.warning('Your job will now be cleaned up.')
        os.remove(jobname + '.pbs')
        self.pbs_log.warning('%s.pbs has been deleted.', jobname)
        os.remove(jobname + '.py')
        self.pbs_log.warning('%s.py has been deleted.' % jobname)

    def format_template_string(self, code=None, template=None, attributes=None):

        if code and attributes is not None:
            code_template = string.Template(code)
            code = code_template.substitute(attributes)
        elif template and attributes is not None:
            with open(template, 'r') as tem:
                code = tem.read()
                code_template = string.Template(code)
                code = code_template.substitute(attributes)
        else:
            code = None
        return code

    def write_template_string(self, code, extension):
        filename = Path(self.pbs_work_dir) / self.job_name + extension
        with open(filename, 'w') as f:
            f.write(code)

    def create_header_section(self, file):
        with open(file, 'a') as f:
            f.write("# Author:  %s\n" % self.author)
            f.write("# Date Created: %s\n" % self.current_date)
            f.write("# Project Name: %s\n" % self.project_name)
            f.write("# Description: %s\n\n" % self.description)

    def create_directives_section(self, file):
        with open(file, 'a') as f:
            f.write("#PBS -l %s\n" % self.resource_str)
            f.write("#PBS -l cput=%s\n" % self.cputime)
            f.write("#PBS -l walltime=%s\n" % self.walltime)
            f.write("#PBS -N %s\n" % self.job_name)
            f.write("#PBS -o %s\n" % (self.job_name + ".o"))
            f.write("#PBS -e %s\n" % (self.job_name + ".e"))
            f.write("#PBS -M %s\n" % self.email)
            if self.directive_list is not None:
                for directive in self.directive_list:
                    f.write("#PBS %s\n" % directive)
            f.write("\n")

    def create_commands_section(self, file):
        with open(file, 'a') as f:
            f.write("cd %s\n" % self.pbs_work_dir)
            if self.activate_script is not None:
                f.write("source %s\n" % self.activate_script)
            if self.primary_pbs_cmd is not None:
                f.write(self.primary_pbs_cmd + "\n")
            if self.pbs_command_list:
                for cmd in self.pbs_command_list:
                    f.write(cmd + "\n")
            if self.email_command is not None:
                f.write(self.email_command + "\n")

    def create_scripts(self, python_code=None, python_template=None, python_attributes=None,
                       pbs_code=None, pbs_template=None, pbs_attributes=None, custom_pbs=False):
        """
        Python code is supplied as a file or as a string.  The string or file can be a template string or file, which
        can be supplied a dictionary to fill in the templated variables.  Additionally, a PBS file (the default from the
        package or a custom) is used in conjunction with the python file in order to submit the job.  The PBS file can
        be generated from a template or a custom script can be generated with a minimum set of parameters already set
        up.

        :param python_code:
        :type python_code:
        :param python_template:
        :type python_template:
        :param python_attributes:
        :type python_attributes:
        :param pbs_code:
        :type pbs_code:
        :param pbs_template:
        :type pbs_template:
        :param pbs_attributes:
        :type pbs_attributes:
        :param custom_pbs:
        :type custom_pbs:
        :return:
        :rtype:
        """

        # Configure the Python Code
        python_code = self.format_template_string(code=python_code, template=python_template, attributes=python_attributes)
        if python_code is None:
            raise ValueError("The python code string or template file needs to be given.")
        self.write_template_string(python_code, extension=".py")
        # Configure the PBS Code
        if not custom_pbs:
            pbs_code = self.format_template_string(code=pbs_code, template=pbs_template, attributes=pbs_attributes)
            if pbs_code is None:
                pbs_code = self.format_template_string(template=self.pbs_template, attributes=pbs_attributes)
                if pbs_code is None:
                    raise ValueError("Please supply pbs attributes.")
            self.write_template_string(pbs_code, extension=".pbs")
        else:
            self.create_header_section(file=self.pbs_script)
            self.create_directives_section(file=self.pbs_script)
            self.create_commands_section(file=self.pbs_script)

