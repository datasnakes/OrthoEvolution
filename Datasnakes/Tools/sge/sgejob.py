import contextlib
from subprocess import run, CalledProcessError, PIPE
import os
from pkg_resources import resource_filename
import json


from Datasnakes.Tools.logit import LogIt
from Datasnakes.Tools.sge import basejobids, writecodefile, import_temp, file2str
from Datasnakes.Tools.sge.sgeconfig import __DEFAULT__, __CUSTOM__
from Datasnakes.Manager.config import templates


class BaseSGEJob(object):
    """Base class for simple jobs."""
    def __init__(self, length, base_jobname):
        """Initialize job attributes."""
        self.default_job_attributes = __DEFAULT__
        self.custom_job_attributes = __CUSTOM__
        self.file2str = file2str
        self.job_id, self.jobname = self._configure(length, base_jobname)
        self.sgejob_log = LogIt().default(logname="SGE JOB", logfile=None)

        # Import the temp.pbs file using pkg_resources
        self.temp_pbs = resource_filename(templates.__name__, "temp.pbs")

    @classmethod
    def _configure(cls, length, base_jobname):
        """Configure job attributes or set it up."""
        baseid, base = basejobids(length, base_jobname)
        return baseid, base

    @classmethod
    def updatedict(dictionary, new_attributes):
        dict_input = json.dumps(dictionary)
        dict_input.format(new_attributes)
        dict_output = json.loads(dict_input)

        return dict_output

    @classmethod
    def _cleanup(cls, jobname):
        """Clean up job scripts."""
        os.remove(jobname + '.pbs')
        os.remove(jobname + '.py')


class SGEJob(BaseSGEJob):
    """Create a qsub/pbs job & script for the job to execute."""
    def __init__(self, base_jobname):
        super().__init__(length=5, base_jobname=base_jobname)

    def debug(self, code):
        """Debug the SGEJob."""
        self.submit(code, debug=True)

    def submit(self, code, cleanup=False, default=True, debug=False):
        """Create and submit a qsub job.

        Submit python code."""
        if debug:
            self.sgejob_log.info('You decided to debug your SGEJob')

        # TIP If python is in your environment as only 'python' update that.
        # If default, a python file will be created from code that is used.
        # Allow user input to be a python file
        if os.path.isfile(code) and str(code).endswith('.py'):
            code_str = self.file2str(code)
            self.sgejob_log.info('%s converted to string.' % code)
        elif type(code) == str:
            code_str = code

        if default:
            self.sgejob_log.info('You are running a job with default attributes.')
            writecodefile(filename=self.jobname, code=code_str, language='python')
            pyfilename = self.jobname + '.py'
            self.sgejob_log.info('%s python file has been created.' % pyfilename)

            # Create the pbs script from the template or dict
            pbstemp = import_temp(self.temp_pbs)

            attributes = self.default_job_attributes
            pbsfilename = self.jobname + '.pbs'
            with open(pbsfilename, 'w') as pbsfile:
                pbsfile.write(pbstemp.substitute(attributes))
                pbsfile.close()
            self.sgejob_log.info('%s has been created.' % pbsfilename)
        else:
            msg = 'Custom SGEJob creation is not yet implemented.'
            self.sgejob_log.error(msg)
            raise NotImplementedError(msg)
            # TODO Add custom job creation

        if debug:
            with contextlib.suppress(CalledProcessError):
                cmd = 'qsub ' + self.jobname + '.pbs'  # this is the command
                self.sgejob_log.info('`%s` will be run.' % cmd)
                self.sgejob_log.info('Job submitted.')
                self.cleanup(self.jobname)

        else:
            with contextlib.suppress(CalledProcessError):
                cmd = ['qsub ' + self.jobname + '.pbs']  # this is the command
                # Shell MUST be True
                cmd_status = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

                if cmd_status.returncode == 0:  # Command was successful.
                    self.sgejob_log.info('Job submitted.')
                    # TODO add a check to for job errors or check for error file.

                else:  # Unsuccessful. Stdout will be '1'
                    self.sgejob_log.error('PBS job not submitted.')
        # TODO Introduce a wait function w/qstat to check for job completion
