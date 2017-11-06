import contextlib
from subprocess import run, CalledProcessError, PIPE
import os

from Datasnakes.Tools.sge import basejobids, writecodefile, import_temp, file2str
from Datasnakes.Tools.sge.sgeconfig import __DEFAULT__, __CUSTOM__


class BaseSGEJob(object):
    """Base class for simple jobs."""
    def __init__(self, jobname):
        """Initialize job attributes."""
        self.default_job_attributes = __DEFAULT__
        self.custom_job_attributes = __CUSTOM__
        self.file2str = file2str
        self.job_id, self.base_jobname = self._configure(jobname)

    def _configure(self, jobname):
        """Configure job attributes or set it up."""
        baseid, base = basejobids(jobname)
        return baseid, base

    @classmethod
    def _cleanup(jobname):
        """Clean up job scripts."""
        os.remove(jobname + '.pbs')
        os.remove(jobname + '.py')

    def dryrun(self):
        pass


class SGEJob(BaseSGEJob):
    """Create multiple jobs & scripts for each job to run."""
    def __init__(self, jobname):
        super().__init__(jobname)

    def submit(self, code, cleanup=False, default=True):
        """Create and submit a qsub job.

        Submit python code."""
        # TIP If python is in your environment as only 'python' update that.
        # TODO-SDH add a default parameter or custom parameter
        # If default, a python file will be created from code that is used.

        # Allow user input to be a python file
        if os.path.isfile(code) and str(code).endswith('.py'):
            code_str = self.file2str(code)
        elif type(code) == str:
            code_str = code

        if default:
            writecodefile(filename=self.base_jobname, code=code_str,
                          language='python')

            # Create the pbs script from the template or dict
            pbstemp = import_temp('temp.pbs')

            attributes = self.default_job_attributes

            with open(self.base_jobname + '.pbs', 'w') as pbsfile:
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
            cmd = ['qsub ' + self.base_jobname + '.pbs']  # this is the command
            # Shell MUST be True
            cmd_status = run(cmd, stdout=PIPE, stderr=PIPE, shell=True)

            if cmd_status.returncode == 0:  # Command was successful.
                print('Job submitted.\n')
                # TODO add a check to for job errors or check for error file.

            else:  # Unsuccessful. Stdout will be '1'
                print("PBS job not submitted.")

        # TODO Introduce a wait function w/qstat to check for job completion
