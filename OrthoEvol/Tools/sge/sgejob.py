from subprocess import run, CalledProcessError, PIPE
import os
from time import sleep
from pkg_resources import resource_filename

from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Tools.sge import basejobids, writecodefile, import_temp, file2str
from OrthoEvol.Tools.sge.sgeconfig import __DEFAULT__
from OrthoEvol.Manager.config import templates
from OrthoEvol.Tools.sge import Qstat


class BaseSGEJob(object):
    """Base class for simple jobs."""
    def __init__(self, base_jobname):
        """Initialize job attributes."""
        self.base_jobname = base_jobname
        self.default_job_attributes = __DEFAULT__
        self.file2str = file2str
        self.sgejob_log = LogIt().default(logname="SGE JOB", logfile=None)
        self.pbsworkdir = os.getcwd()

        # Import the temp.pbs file using pkg_resources
        self.temp_pbs = resource_filename(templates.__name__, "temp.pbs")

    @classmethod
    def _configure(cls, length, base_jobname):
        """Configure job attributes or set it up."""
        baseid, base = basejobids(length, base_jobname)
        return baseid, base

    def debug(self, code):
        """Debug the SGEJob."""
        raise NotImplementedError

    def _cleanup(self, jobname):
        """Clean up job scripts."""
        self.sgejob_log.warning('Your job will now be cleaned up.')
        os.remove(jobname + '.pbs')
        self.sgejob_log.warning('%s.pbs has been deleted.', jobname)
        os.remove(jobname + '.py')
        self.sgejob_log.warning('%s.py has been deleted.' % jobname)

    def wait_on_job_completion(self, job_id):
        """Use Qstat to monitor your job."""
        # TODO Allow either slack notifications or email or text.
        qwatch = Qstat().watch(job_id)
        if qwatch == 'Job id not found.':
            self.sgejob_log.info('%s has finished.' % job_id)
            sleep(30)
        elif qwatch == 'Waiting for %s to start running.' % job_id:
            self.sgejob_log.info('%s is queued to run.' % job_id)
            self.sgejob_log.info('Waiting for %s to start.' % job_id)
            sleep(30)
            self.wait_on_job_completion(job_id)
        elif qwatch == 'Waiting for %s to finish running.' % job_id:
            self.sgejob_log.info('%s is running.' % job_id)
            self.sgejob_log.info('Waiting for %s to finish.' % job_id)
            sleep(30)
            self.wait_on_job_completion(job_id)
        else:
            self.wait_on_job_completion(job_id)

    def submitjob(self, cleanup, wait=True):
        """Submit a job using qsub."""
        try:
            cmd = ['qsub ' + self.jobname + '.pbs']  # this is the command
            # Shell MUST be True
            cmd_status = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, check=True)
        except CalledProcessError as err:
            self.sgejob_log.error(err.stderr.decode('utf-8'))
            if cleanup:
                self._cleanup(self.jobname)
        else:
            if cmd_status.returncode == 0:  # Command was successful.
                # The cmd_status has stdout that must be decoded.
                # When a qsub job is submitted, the stdout is the job id.
                submitted_jobid = cmd_status.stdout.decode('utf-8')
                self.sgejob_log.info(self.jobname + ' was submitted.')
                self.sgejob_log.info('Your job id is: %s' % submitted_jobid)
                if wait is True:
                    self.wait_on_job_completion(submitted_jobid)
                self._cleanup(self.jobname)

            else:  # Unsuccessful. Stdout will be '1'
                self.sgejob_log.error('PBS job not submitted.')


class SGEJob(BaseSGEJob):
    """Create a qsub/pbs job & script to submit python code."""
    def __init__(self, email_address, base_jobname=None):
        super().__init__(base_jobname=base_jobname)
        self.email = email_address
        self.attributes = self.default_job_attributes
        self.jobname = self.default_job_attributes['job_name']
        if base_jobname is not None:
            _, self.jobname = self._configure(base_jobname=base_jobname,
                                              length=5)
            self.attributes = self._update_default_attributes()

    def _update_default_attributes(self):
        pyfile_path = os.path.join(self.pbsworkdir, self.jobname + '.py')
        # These attributes are automatically updated if a jobname is given.
        new_attributes = {'email': self.email,
                          'job_name': self.jobname,
                          'outfile': self.jobname + '.o',
                          'errfile': self.jobname + '.e',
                          'script': self.jobname,
                          'log_name': self.jobname + '.log',
                          'cmd': 'python3.6 ' + pyfile_path,
                          }
        self.default_job_attributes.update(new_attributes)

        return self.default_job_attributes

    def submit_pycode(self, code, wait=True, cleanup=True, default=True):
        """Create and submit a qsub job.

        Submit python code."""
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
            self.sgejob_log.info('%s has been created.' % pyfilename)

            # Create the pbs script from the template or dict
            pbstemp = import_temp(self.temp_pbs)
            pbsfilename = self.jobname + '.pbs'

            with open(pbsfilename, 'w') as pbsfile:
                pbsfile.write(pbstemp.substitute(self.attributes))
                pbsfile.close()
            self.sgejob_log.info('%s has been created.' % pbsfilename)

            # Submit the job using qsub
            self.submitjob(cleanup=cleanup, wait=wait)
        else:
            msg = 'Custom SGEJob creation is not yet implemented.'
            raise NotImplementedError(msg)
            # TODO Add custom job creation

        # Submit the job using qsub
        self.submitjob(cleanup=cleanup, wait=wait)


class MultiSGEJobs(SGEJob):
    """Create multiple qsub/pbs jobs."""
    def __init__(self, email_address):
        super().__init__(email_address=email_address, base_jobname='multi')
        _, self.jobname = self._configure(base_jobname=self.base_jobname,
                                          length=5)
        self.attributes = self._update_default_attributes()
