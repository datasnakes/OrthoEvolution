"""Access qstat information about SGE jobs."""
from subprocess import check_output, CalledProcessError
import getpass
import re

class Qstat(object):
    def __init__(self):
        """Initialize class."""
        _username = getpass.getuser()
        self.username = _username
        self.split_regex = re.compile(r'\s+')
        
    def qstatinfo(self, qstat_path='qstat'):
        """Retrieve qstat output."""
        try:
            qstatinfo = check_output([qstat_path])
        except CalledProcessError as cpe:
            return_code = 'qstat returncode: %s' % cpe.returncode
            std_error = 'qstat standard output: %s' % cpe.stderr
            print(return_code + '\n' + std_error)
        except FileNotFoundError:
            raise FileNotFoundError('qstat is not on your machine.')

        jobs = self._output_parser(qstatinfo)
        
        return jobs

#        self.running_jobs = sum(j['status'] in ['R', 'Q'] for j in jobs)
#
#        return self.running_jobs

# TODO Create function for getting a list or dict of running jobs.
# TODO Create function for getting a list or dict of queued jobs.
# TODO Create function for getting a list or dict of running jobs for a user.
# TODO Create a functions that checks every few minutes until job finishes.

    def _output_parser(self, output):
        """Parse output from qstat pbs commandline program."""
        lines = output.decode('utf-8').split('\n')
        del lines[:5]
        jobs = []
        for line in lines:
            els = self.split_regex.split(line)
            try:
            	j = {"job_id": els[0], "name": els[1], "user": els[2], "elapsed_time": els[3],
                 	"status": els[4], "queue": els[5]}
            	jobs.append(j)

            except IndexError:
                pass

        return jobs
