"""Access a list of SGE jobs."""
from subprocess import run, CalledProcessError, PIPE
import getpass


class Qstat(object):
    def __init__(self):
        """Initialize class."""
        _username = getpass.getuser()
        self.username = _username
        self.qstatinfo()


    def qstatinfo(self, qstat_path='qstat', option='-u'):
        """Retrieve qstat output."""
        try:
            qstatinfo = run([qstat_path, option, self.username],
                            stdout=PIPE, stderr=PIPE, shell=True)
        except CalledProcessError as cpe:
            return_code = 'qstat returncode: %s' % cpe.returncode
            std_error = 'qstat standard output: %s' % cpe.stderr
            print(return_code + '\n' + std_error)
        except FileNotFoundError:
            raise FileNotFoundError('qstat is not on your machine.')

        jobs = self._output_parser(qstatinfo.stdout)

        self.running_jobs = sum(j['status'] in ['R', 'Q'] for j in jobs)

        return self.running_jobs

# TODO Add function that will parse output of qstat
    def _output_parser(self, output):
        """
        Parse output from qstat pbs commandline program

        Arguments:
        ----------
        output: str, output obtained from qstat command
        Returns:
        ---------
        list, of all jobs. Each job is represented as a dictionary conainting
            all relevant information.
        """
        lines = output.split('\n')
        del lines[:5]
        jobs = []
        for line in lines:
            els = self.split_regex.split(line)
            try:
            	j = {"id_": els[0], "user": els[1], "queue": els[2], "name": els[3],
                 	"status": els[9], "elapsed_time": els[10]}
            	jobs.append(j)

            except IndexError:
                pass

        return jobs
