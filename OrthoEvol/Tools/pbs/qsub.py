import getpass
import string
import random
import os
import shutil
import subprocess as sp
from pathlib import Path
from collections import OrderedDict
from pkg_resources import resource_filename
from datetime import datetime as d
from time import sleep
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Manager.config import templates
from OrthoEvol.utilities import FullUtilities


class BaseQsub(object):
    """Base class for PBS jobs."""

    def __init__(self, pbs_script, job_name, base_job_id=None, pbs_working_dir=None):

        self.qsub_log = LogIt().default(logname="PBS - QSUB", logfile=None)
        self.qsub_utils = FullUtilities()
        self.base_name = job_name
        if base_job_id is None:
            self.base_id, self.job_name = self.get_base_job_name()
        else:
            self.base_id = base_job_id
            self.job_name = self.base_name + "_%s" % base_job_id
        # PBS working directory
        if pbs_working_dir is None:
            self.pbs_working_dir = Path(os.getcwd()) / Path(self.job_name)
        else:
            self.pbs_working_dir = Path(pbs_working_dir) / Path(self.job_name)
        # PBS script
        self.pbs_script = Path(self.pbs_working_dir) / (self.job_name + '.pbs')
        if pbs_script is not None:
            self.supplied_pbs_script = Path(pbs_script)
        else:
            self.supplied_pbs_script = self.pbs_script

        if not self.pbs_working_dir.exists():
            self.pbs_working_dir.mkdir(parents=True)

        self.pbs_job_id = None
        self.qsub_job_directory = None

    def get_base_job_name(self, length=5):
        base_id = ''.join(random.sample(string.ascii_letters + string.digits, length))
        job_name = self.base_name + "_%s" % base_id

        return base_id, job_name

    def copy_supplied_script(self, supplied_script, new_script):
        # if the supplied script exists make sure its in the PBS working directory
        if Path(supplied_script).exists():
            if not Path(supplied_script) == Path(new_script):
                shutil.copy(str(supplied_script), str(new_script))
        else:
            raise FileExistsError("The PBS script does not exists.")

    def submit_pbs_script(self, cmd=None):
        """Submit a job using qsub.

        :param cleanup: (Default value = False)
        :param wait: (Default value = True)
        """

        try:
            if cmd is None:
                cmd = ['qsub ' + str(self.pbs_script)]  # this is the command
            # Shell MUST be True
            proc = self.qsub_utils.system_cmd(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        except sp.CalledProcessError as err:
            self.qsub_log.error(err.stderr)
        else:
            if proc.returncode == 0:  # Command was successful.
                # The cmd_status has stdout that must be decoded.
                # When a qsub job is submitted, the stdout is the job id.
                self.pbs_proc = proc
                self.pbs_job_id = proc.stdout[0]
                self.qsub_log.info(str(self.pbs_script) + ' was submitted.')
                self.qsub_log.info('Your job id is: %s' % self.pbs_job_id)
                self.qsub_job_directory = self.pbs_working_dir / self.pbs_job_id

            else:  # Unsuccessful. Stdout will be '1'
                self.qsub_log.error('PBS job not submitted.')


class Qsub(BaseQsub):

    def __init__(self, python_script=None, author=None, project_name="OrthoEvol", description="This is a basic pbs job",
                 date_format='%a %b %d %I:%M:%S %p %Y', chunk_resources=None, cput='72:00:00', walltime='48:00:00',
                 email=None, directive_list=None, pbs_command_list=None,
                 email_command=None, **kwargs):

        super().__init__(**kwargs)

        # Python script
        self.python_script = Path(self.pbs_working_dir) / (self.job_name + '.py')
        if python_script is not None:
            self.supplied_python_script = Path(python_script)
        else:
            self.supplied_python_script = self.python_script

        self.pbs_template = resource_filename(templates.__name__, "temp.pbs")

        # PBS - header info
        if author is None:
            self.author = getpass.getuser()
        else:
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
        self.pbs_command_list = []
        self.pbs_command_list.extend(pbs_command_list)
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
                "mem": "6gb"
            })

        resource_list = []
        for k, v in chunk_resources.items():
            if v is not None:
                resource_list.append("%s=%s" % (k, v))
        resource_str = ":".join(resource_list)
        return resource_str

    def format_template_string(self, code=None, template=None, attributes=None):

        if code and attributes is not None:
            code_template = string.Template(code)
            code = code_template.substitute(attributes)
        elif template and attributes is not None:
            with open(str(template), 'r') as tem:
                code = tem.read()
                code_template = string.Template(code)
                code = code_template.substitute(attributes)
        else:
            code = None
        return code

    def write_template_string(self, code, extension, file=None):
        if file is None:
            filename = Path(self.pbs_working_dir) / str(self.job_name + extension)
        else:
            filename = Path(file)
        with open(str(filename), 'w') as f:
            f.write(code)

    def create_header_section(self, file):
        with open(file, 'a') as f:
            f.write("# Author:  %s\n" % self.author)
            f.write("# Date Created: %s\n" % self.current_date)
            f.write("# Project Name: %s\n" % self.project_name)
            f.write("# Description: %s\n\n" % self.description)

    def create_directives_section(self, file):
        with open(file, 'a') as f:
            if self.resource_str is not None:
                f.write("#PBS -l %s\n" % self.resource_str)
            if self.cputime is not None:
                f.write("#PBS -l cput=%s\n" % self.cputime)
            if self.walltime is not None:
                f.write("#PBS -l walltime=%s\n" % self.walltime)
            if self.job_name is not None:
                f.write("#PBS -N %s\n" % self.job_name)
                f.write("#PBS -o %s\n" % str(self.pbs_working_dir / (self.job_name + ".o")))
                f.write("#PBS -e %s\n" % str(self.pbs_working_dir / (self.job_name + ".e")))
            if self.email is not None:
                f.write("#PBS -M %s\n" % self.email)
            if self.directive_list is not None:
                for directive in self.directive_list:
                    f.write("#PBS %s\n" % directive)
            f.write("\n")

    def create_commands_section(self, file):
        with open(file, 'a') as f:
            f.write("cd %s\n" % str(self.pbs_working_dir))
            if self.pbs_command_list:
                for cmd in self.pbs_command_list:
                    f.write(cmd + "\n")
            if self.email_command is not None:
                f.write(self.email_command + "\n")

    def format_python_script(self, py_template_string=None, py_template_file=None, python_attributes=None):
        # Configure the Python Code
        python_code = self.format_template_string(code=py_template_string, template=py_template_file,
                                                  attributes=python_attributes)
        if python_code is not None:
            self.write_template_string(python_code, file=self.supplied_python_script)

    def set_up_pbs_script(self, pbs_template_string=None, pbs_template_file=None, pbs_attributes=None):
        """
        Python code is supplied as a file or as a string.  The string or file can be a template string or file, which
        can be supplied a dictionary to fill in the templated variables.  Additionally, a PBS file (the default from the
        package or a custom) is used in conjunction with the python file in order to submit the job.  The PBS file can
        be generated from a template or a custom script can be generated with a minimum set of parameters already set
        up.

        :param pbs_template_string:
        :type pbs_template_string:
        :param pbs_template_file:
        :type pbs_template_file:
        :param pbs_attributes:
        :type pbs_attributes:
        :return:
        :rtype:
        """

        # Configure the PBS Code
        if pbs_attributes is not None:
            if pbs_template_file is not None:
                pbs_code = self.format_template_string(template=pbs_template_file, attributes=pbs_attributes)
            elif pbs_template_string is not None:
                pbs_code = self.format_template_string(code=pbs_template_string, attributes=pbs_attributes)
            else:
                raise ValueError("Please supply the pbs_template_file or pbs_template_string to generate the proper"
                                 "pbs script.")
            self.write_template_string(pbs_code, file=self.supplied_pbs_script)
        else:
            self.create_header_section(file=self.supplied_pbs_script)
            self.create_directives_section(file=self.supplied_pbs_script)
            self.create_commands_section(file=self.supplied_pbs_script)

    def submit_python_job(self, cmd=None, py_template_string=None, py_template_file=None, python_attributes=None,
                          pbs_template_string=None, pbs_template_file=None, pbs_attributes=None,
                          custom_python_cmd=None, rerun=False):
        if not rerun:
            self.format_python_script(py_template_string=py_template_string, py_template_file=py_template_file,
                                      python_attributes=python_attributes)
            if custom_python_cmd is not None:
                self.pbs_command_list.append(custom_python_cmd)
            else:
                self.pbs_command_list.append("python %s" % self.python_script)
            self.copy_supplied_script(supplied_script=self.supplied_python_script, new_script=self.python_script)
            self.set_up_pbs_script(pbs_template_string=pbs_template_string, pbs_template_file=pbs_template_file,
                                   pbs_attributes=pbs_attributes)
            self.copy_supplied_script(supplied_script=self.supplied_pbs_script, new_script=self.pbs_script)
        else:
            self.submit_pbs_script(cmd=cmd)
