"""Access a list of SGE jobs."""
from subprocess import run, CalledProcessError, PIPE
import getpass


class Qstat(object):
    def __init__(self):
        """
        Parameters
        ----------
        qstat_path : string
            The path to the qstat executable.
        Returns
        -------
        queue_info : list
            A list of jobs in 'queue_info'. Jobs are dictionaries with both
            string keys and string names.
        job_info : list
            A list of jobs in 'job_info'.
        """
        _username = getpass. getuser()
        self.username = _username
        self.qstatinfo()
        

    @classmethod
    def qstatinfo(cls, qstat_path='qstat', option=''):
        """
        Parameters
        ----------
        qstat_path : string
            The path to the qstat executable.
        Returns
        -------
        qstatxml : string
            The xml stdout string of the 'qstat -xml' call.
        """
        try:
            qstatinfo = run([qstat_path, option], stdout=PIPE, stderr=PIPE, shell=True)
        except CalledProcessError as cpe:
            return_code = 'qstat returncode: %s' % cpe.returncode
            std_error = 'qstat standard output: %s' % cpe.stderr
            print(return_code + '\n' + std_error)
        except FileNotFoundError:
            raise FileNotFoundError('qstat is not on your machine.')
        
        return qstatinfo.stdout

# TODO Add function that will parse output of qstat