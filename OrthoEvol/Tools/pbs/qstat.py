"""Access qstat information about SGE jobs."""
from subprocess import check_output, CalledProcessError
import getpass
import re
from time import sleep

class Qstat(object):
    """Access and parse information from qstat commands."""

    def __init__(self):
        """Get user name."""

        _username = getpass.getuser()
        self.username = _username
        self.split_regex = re.compile(r'\s+')

    def qstatinfo(self, qstat_cmd='qstat'):
        """Retrieve qstat output.

        :param qstat_cmd: Qstat command to run.  (Default value = 'qstat')
        """

        try:
            qstatinfo = check_output([qstat_cmd])
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

        :param output: Input the output of the qstat command.
        :return: A list of dictionaries for each job.
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
        """Wait until a job or list of jobs finishes and get updates.

        :param job_id: Input a job id to be monitored.
        :return: Status of job id.
        """

        jobs = self.qstatinfo()
        rids = [j['job_id'] for j in jobs if j['user'] == self.username
                and j['status'] == 'R']
        qids = [j['job_id'] for j in jobs if j['user'] == self.username
                and j['status'] == 'Q']
        if job_id in qids:
            yield 'Waiting for %s to start running.' % job_id
            self.watch(job_id)
        elif job_id in rids:
            yield 'Waiting for %s to finish running.' % job_id
            self.watch(job_id)
        else:
            return 'Job id not found.'

    def wait_on_job_completion(self, job_id):
        """Use Qstat to monitor your job.

        :param job_id: The job id to be monitored.
        """

        # TODO Allow either slack notifications or email or text.
        qwatch = self.watch(job_id)
        if qwatch == 'Job id not found.':
            #self.pbs_log.info('%s has finished.' % job_id)
            sleep(30)
        elif qwatch == 'Waiting for %s to start running.' % job_id:
            #self.pbs_log.info('%s is queued to run.' % job_id)
            #self.pbs_log.info('Waiting for %s to start.' % job_id)
            sleep(30)
            self.wait_on_job_completion(job_id)
        elif qwatch == 'Waiting for %s to finish running.' % job_id:
            #self.pbs_log.info('%s is running.' % job_id)
            #self.pbs_log.info('Waiting for %s to finish.' % job_id)
            sleep(30)
            self.wait_on_job_completion(job_id)
        else:
            self.wait_on_job_completion(job_id)