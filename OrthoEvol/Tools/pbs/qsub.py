import getpass
import string
import random
import os
import subprocess as sp
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
            cmd = ['qsub ' + self.pbs_script]  # this is the command
            # Shell MUST be True
            cmd_status = self.pbs_utils.system_cmd(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, check=True)
        except sp.CalledProcessError as err:
            self.pbs_log.error(err.stderr.decode('utf-8'))
        else:
            if cmd_status.returncode == 0:  # Command was successful.
                # The cmd_status has stdout that must be decoded.
                # When a qsub job is submitted, the stdout is the job id.
                submitted_jobid = cmd_status.stdout.decode('utf-8')
                self.pbs_log.info(self.pbs_script + ' was submitted.')
                self.pbs_log.info('Your job id is: %s' % submitted_jobid)

            else:  # Unsuccessful. Stdout will be '1'
                self.pbs_log.error('PBS job not submitted.')


class Qsub(BaseQsub):

    def __init__(self, pbs_script=None, author=getpass.getuser(), project_name="OrthoEvol", description="This is a basic pbs job",
                 date_format='%a %b %d %I:%M:%S %p %Y', chunk_resources=None, cput='72:00:00', walltime='48:00:00',
                 job_name=None, pbs_work_dir=None, script_cmd=None, email=None, directive_list=None,
                 activate_script=None):

        super().__init__(pbs_script=pbs_script)

        self.pbs_template = resource_filename(templates.__name__, "temp.pbs")
        # Get path to activate script for python virtual environment
        self.activate_script = activate_script
        # Set up commented script header
        self.author = author
        self.project_name = project_name
        self.description = description
        self.current_date = d.now().strftime(date_format)

        # Set up PBS directives/attributes
        # resources
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
        self.resource_str = ":".join(resource_list)
        self.cputime = cput
        self.walltime = walltime
        # other attributes
        self.base_name = job_name
        self.base_id, self.job_name = self.get_base_job_name()
        self.directive_list = directive_list

        # Set up PBS variables
        self.pbs_work_dir = pbs_work_dir
        self.script_cmd = script_cmd
        self.email = email

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

    def write_template_string(self):
        pass

    def submit_python_script(self, python_code=None, python_template=None, python_attributes=None,
                             pbs_code=None, pbs_template=None, pbs_attributes=None):
        """
        Python code is supplied as a file or as a string.  The string or file can be a template string or file, which
        can be supplied a dictionary to fill in the templated variables.  Additionally, a PBS file (the default from the
        package or a custom) is used in conjunction with the python file in order to submit the job.

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
        :return:
        :rtype:
        """

        # Configure the Python Code
        python_code = self.format_template(code=python_code, template=python_template, attributes=python_attributes)
        if python_code is None:
            raise ValueError("The python code string or template file needs to be given.")
        pbs_code = self.format_template(code=pbs_code, template=pbs_template, attributes=pbs_attributes)
        if pbs_code is None:
            pbs_code = self.format_template(template=self.pbs_template, attributes=pbs_attributes)
            if pbs_code is None:
                raise ValueError("Please supply pbs attributes.")
