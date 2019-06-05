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
        self.temp_pbs = resource_filename(templates.__name__, "temp.pbs")
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
