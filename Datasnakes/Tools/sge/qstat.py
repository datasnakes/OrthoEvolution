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

    def _output_parser(self, output):
        """Parse output from qstat pbs commandline program.

        Returns a list of dictionaries for each job.
        """
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

    def all_job_ids(self):
        """Retrieve a list of all jobs running or queued."""
        jobs = self.qstatinfo()
        ids = [j['job_id'] for j in jobs]
        return ids

    def all_running_jobs(self):
        """Retrieve a list of running jobs."""
        jobs = self.qstatinfo()
        ids = [j['job_id'] for j in jobs if j['status'] == 'R']
        return ids

    def all_queued_jobs(self):
        """Retrieve a list of queued jobs."""
        jobs = self.qstatinfo()
        ids = [j['job_id'] for j in jobs if j['status'] == 'Q']
        return ids

    def myjobs(self):
        """Retrieve a list of all the current user's jobs."""
        jobs = self.qstatinfo()
        ids = [j['job_id'] for j in jobs if j['user'] == self.username]
        if len(ids) < 1:
            return 'You have no jobs running or queued.'
        else:
            rids = [j['job_id'] for j in jobs if j['user'] == self.username
                    and j['status'] == 'R']
            qids = [j['job_id'] for j in jobs if j['user'] == self.username
                    and j['status'] == 'Q']
            return 'Running jobs: %s\nQueued jobs: %s' % (rids, qids)

    def watch(self, job_id):
        """Wait until a job or list of jobs finishes and get updates."""
        jobs = self.qstatinfo()
        rids = [j['job_id'] for j in jobs if j['user'] == self.username
                and j['status'] == 'R']
        qids = [j['job_id'] for j in jobs if j['user'] == self.username
                and j['status'] == 'Q']
        if job_id in qids:
            return 'Waiting for %s to start running.' % job_id
            self.watch(job_id)
        elif job_id in rids:
            return 'Waiting for %s to finish running.' % job_id
            self.watch(job_id)
        else:
            return 'Job id not found.'
