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

    def __init__(self, job_name, pbs_script=None, base_job_id=None, pbs_working_dir=None):
        """
        The BaseQsub class is the most basic means to create PBS jobs using qsub.
        It takes a given pbs script and creates a new directory for the pbs job using
        a randomly generating string of letters/numbers.  After the supplied pbs script
        is copied into the new folder it can be submitted to the system.  After submission
        the PBS job id will be available as a class variable.

        :param pbs_script:   The path to the PBS script that will be submitted to the system.
        :type pbs_script:  str.
        :param job_name:  A job name that will be used for naming the pbs working directory, and
        for the -N parameter in the PBS script.
        :type job_name:  str.
        :param base_job_id:  While usually none, this can be used to re-run a PBS job that was
        created by this class.
        :type base_job_id:  str.
        :param pbs_working_dir:  The path to the PBS working directory where the scripts will be
        copied to and where the job will be started and run by the PBS system.
        :type pbs_working_dir:  str.
        """

        self.qsub_log = LogIt().default(logname="PBS - QSUB", logfile=None)
        self.qsub_utils = FullUtilities()
        self.base_name = job_name
        if not base_job_id:
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
        self.pbs_proc = None
        self.qsub_job_directory = None

    def get_base_job_name(self, length=5):
        """
        Returns a tuple of the base job id and the job_name.  The base
        job id can be of any length.

        :param length:  The number of characters used in the base job id.
        :type length:  int.
        """
        base_id = ''.join(random.sample(string.ascii_letters + string.digits, length))
        job_name = self.base_name + "_%s" % base_id

        return base_id, job_name

    def copy_supplied_script(self, supplied_script, new_script):
        """
        Copies a supplied script/file using the new script/file name and
        location unless the supplied script doesn't exist.

        :param supplied_script:  The script to be copied.
        :type supplied_script:  str.
        :param new_script:  The path (<path>/<new-script-name>) to copy the supplied
        script to.
        :type new_script:  str.
        """
        # if the supplied script exists make sure its in the PBS working directory
        if Path(supplied_script).exists():
            if not Path(supplied_script) == Path(new_script):
                shutil.copy(str(supplied_script), str(new_script))
        else:
            raise FileExistsError("The script does not exists.")

    def submit_pbs_script(self, cmd=None, **kwargs):
        """
        A function to submit the PBS job to the system.

        :param cmd:  A command list/string used qsub system call.
        :type cmd:  list or str.
        :param kwargs:  Keyword arguments used in subprocess.Popen.
        :type kwargs:  dict.
        """

        try:
            if cmd is None:
                cmd = ['qsub', str(self.pbs_script)]  # this is the command
            # Shell MUST be True
            proc = self.qsub_utils.system_cmd(cmd, stdout=sp.PIPE, stderr=sp.PIPE, **kwargs)
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
        """
        The Qsub class is the primary means of creating a job.  It adds the additional
        functionality of creating custom PBS jobs by using template PBS scripts or by
        utilizing the it's parameters to create each section of the PBS script (header,
        directives, and commands sections).  It also gives the ability to format/template/submit
        python based PBS jobs.

        :param python_script:  The path to the python script that will be submitted.
        :type python_script:  str.
        :param author:  The name of the person submitting the PBS job.  Used in the header section
        of the PBS script.
        :type author:  str.
        :param project_name:  The project name associated with the PBS job.  Used in the header
        section of the PBS script.
        :type project_name:  str.
        :param description:  A description of the job used in the header section of the PBS script.
        :type description:  str.
        :param date_format:  The date format used in the header section of the PBS script.
        :type date_format:  str.
        :param chunk_resources:  A dictionary of chunk resources used by the -l parameter.  Used in
        the directives section of the PBS script.
        :type chunk_resources:  dict.
        :param cput:  The cpu time used in the directives section of the PBS script.
        :type cput:  str.
        :param walltime:  The walltime used in the directives section of the PBS script.
        :type walltime:  str.
        :param email:  The email of the author or user.  Used in the directives section of the PBS
        script.
        :type email:  str.
        :param directive_list:  A list of other directives to use in the PBS script.
        :type directive_list:  list.
        :param pbs_command_list:  A list of extra PBS commands to add to the PBS script.
        :type pbs_command_list:  list.
        :param email_command:  A custom email command used in the command section of the PBS script.
        :type email_command:  str.
        :param kwargs:  Keyword arguments given to the BaseQsub class.
        :type kwargs:  dict.
        """

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
        """
        Creates a list of chunk resources that will be used in the
        directives section of the custom PBS script.

        :param chunk_resources:  A dictionary of values used in the PBS directive section.
        The default will look like this after it is added to the PBS script:
          '#PBS -l select=1:ncpus=1:mem=6gb'
        :type chunk_resources:  dict.
        :return:  A string that will added to the directive list.
        :rtype:  str.
        """
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
        """
        Using string Templating, this function takes a dictionary of attributes
        that are used to format the given code or template file.
        See https://docs.python.org/3.4/library/string.html#string.Template.

        :param code:  A string of code that is a template.
        :type code:  str.
        :param template:  The path to a file that is a template.
        :type template:  str.
        :param attributes:  A dictionary of attributes for templating the code/file.
        :type attributes:  dict.
        :return:  The template code with the attribute variables filled in.
        :rtype:  str.
        """

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
        """
        This function writes a code string to a file located in the
        PBS working directory.  It uses the job_name and the extension to
        create the file:  <pbs_working_dir>/<job_name>.<extension>

        :param code:  The code that will be written to the file.
        :type code:   str.
        :param extension:  The extension of the file.  Usually '.pbs' or '.py'.
        :type extension:  str.
        :param file:  A path to the file that will be written.
        :type file:  str.
        """

        if not file:
            filename = Path(self.pbs_working_dir) / str(self.job_name + extension)
        else:
            filename = Path(file)
        with open(str(filename), 'w') as f:
            f.write(code)

    def create_header_section(self, file):
        """
        Create the header section of the supplied PBS file.

        :param file:  A path to the PBS file.
        :type file:  str.
        """

        with open(file, 'a') as f:
            f.write("# Author:  %s\n" % self.author)
            f.write("# Date Created: %s\n" % self.current_date)
            f.write("# Project Name: %s\n" % self.project_name)
            f.write("# Description: %s\n\n" % self.description)

    def create_directives_section(self, file):
        """
        Create the directive section of the PBS file.

        :param file:  A path to the PBS file.
        :type file:  str.
        """
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
        """
        Create the command section of the PBS file.

        :param file:  A path to the PBS file.
        :type file:  str.
        """
        with open(file, 'a') as f:
            f.write("cd %s\n" % str(self.pbs_working_dir))
            if self.pbs_command_list:
                for cmd in self.pbs_command_list:
                    f.write(cmd + "\n")
            if not self.email_command:
                f.write(self.email_command + "\n")

    def format_python_script(self, py_template_string=None, py_template_file=None, python_attributes=None):
        """
        Lightly wrapping the format_template_script, this function formats a python script
        as a string using the supplied attributes.  The code/script is then writen to the
        the pbs working directory using the naming convention:
        <pbs_working_dir>/<job_name>.py

        See https://docs.python.org/3.4/library/string.html#string.Template.


        :param py_template_string:  Python code as a string.
        :type py_template_string:  str.
        :param py_template_file:  A path to the python template.
        :type py_template_file:  str.
        :param python_attributes:  A dictionary of attributes that will be used with the template.
        :type python_attributes:  dict.
        """

        if py_template_file == self.python_script:
            raise FileExistsError("The script provided already exists.  Do not overwrite.")
        # Configure the Python Code
        python_code = self.format_template_string(code=py_template_string, template=py_template_file,
                                                  attributes=python_attributes)

        if python_code is not None:
            self.write_template_string(python_code, file=self.python_script)

    def set_up_pbs_script(self, pbs_template_string=None, pbs_template_file=None, pbs_attributes=None):
        """
        Python code is supplied as a file or as a string.  The string or file can be a template string or file, which
        can be supplied a dictionary to fill in the templated variables.  Additionally, a PBS file (the default from the
        package or a custom) is used in conjunction with the python file in order to submit the job.  The PBS file can
        be generated from a template or a custom script can be generated with a minimum set of parameters already set
        up.

        See https://docs.python.org/3.4/library/string.html#string.Template.

        :param pbs_template_string:  The PBS script as a string.
        :type pbs_template_string:  str.
        :param pbs_template_file:  A path to the PBS template.
        :type pbs_template_file:  str.
        :param pbs_attributes:  A dictionary of attributes that will be used with the template.
        :type pbs_attributes:  dict.
        """

        if pbs_template_file == self.pbs_script:
            raise FileExistsError("The script provided already exists.  Do not overwrite.")

        # Configure the PBS Code
        if pbs_attributes is not None:
            if pbs_template_file is not None:
                pbs_code = self.format_template_string(template=pbs_template_file, attributes=pbs_attributes)
            elif pbs_template_string is not None:
                pbs_code = self.format_template_string(code=pbs_template_string, attributes=pbs_attributes)
            else:
                raise ValueError("Please supply the pbs_template_file or pbs_template_string to generate the proper"
                                 "pbs script.")
            if pbs_code is not None:
                self.write_template_string(pbs_code, file=self.pbs_script)
        else:
            self.create_header_section(file=self.pbs_script)
            self.create_directives_section(file=self.pbs_script)
            self.create_commands_section(file=self.pbs_script)

    def submit_python_job(self, cmd=None, py_template_string=None, py_template_file=None, python_attributes=None,
                          pbs_template_string=None, pbs_template_file=None, pbs_attributes=None,
                          custom_python_cmd=None, rerun=False):
        """
        This function is the primary means of submitting a PBS job that runs a python script.
        It will format and copy the PBS script and the Python script that will be submitted
        by the PBS system.  The PBS script will either be customized using the class parameters
        or it will be formatted using the template parameters.  It also takes custom python commands
        as a list, which is useful for certain actions like activating virutal environments.

        :param cmd:  The command used to run the qsub job.
        :type cmd:  str.
        :param py_template_string:  Python code as a string.
        :type py_template_string:  str.
        :param py_template_file:  A path to the python template.
        :type py_template_file:  str.
        :param python_attributes:  A dictionary of attributes that will be used with the template.
        :type python_attributes:  dict.
        :param pbs_template_string:  The PBS script as a string.
        :type pbs_template_string:  str.
        :param pbs_template_file:  A path to the PBS template.
        :type pbs_template_file:  str.
        :param pbs_attributes:  A dictionary of attributes that will be used with the template.
        :type pbs_attributes:  dict.
        :param custom_python_cmd:  A list of custom commands that are used to create the command
        section of the PBS script.
        :type custom_python_cmd:  list.
        :param rerun:  A flag used to rerun a previously created job.
        :type rerun:  bool.
        """
        if not rerun:
            # Format or copy the python script.
            if python_attributes is not None:
                self.format_python_script(py_template_string=py_template_string, py_template_file=py_template_file,
                                          python_attributes=python_attributes)
            elif not self.python_script.exists():
                self.copy_supplied_script(supplied_script=self.supplied_python_script, new_script=self.python_script)

            # Set up the custom python command
            if custom_python_cmd is not None:
                self.pbs_command_list.append(custom_python_cmd)
            else:
                self.pbs_command_list.append("python %s" % self.python_script)

            # Format the PBS script, create a custom PBS script, or copy the PBS script
            self.set_up_pbs_script(pbs_template_string=pbs_template_string, pbs_template_file=pbs_template_file,
                                   pbs_attributes=pbs_attributes)
            if not self.pbs_script.exists():
                self.copy_supplied_script(supplied_script=self.supplied_pbs_script, new_script=self.pbs_script)

            # Submit job
            self.submit_pbs_script(cmd=cmd)
        else:
            # Submit job as rerun
            self.submit_pbs_script(cmd=cmd)
