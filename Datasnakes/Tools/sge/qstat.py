"""Access a list of SGE jobs."""

from subprocess import check_output, STDOUT, CalledProcessError
import xmltodict


class Qstat(object):
    def __init__(self, xml_option='xml'):
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
        xml = self._qstat2xml(qstat_path='qstat', xml_option='-' + xml_option)
        self.xml = xml
        queue_info, job_info = self._qstatinfo(self.xml)

    def _qstat2xml(self, qstat_path='qstat', xml_option='-xml'):
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
            qstatxml = check_output([qstat_path, xml_option], stderr=STDOUT)
        except CalledProcessError as cpe:
            print('qstat returncode: ', cpe.returncode)
            print('qstat standard output: ', cpe.output)
            raise
        except FileNotFoundError as fnfe:
            fnfe.message = 'Maybe qstat -xml is not on your machine.'
            raise
        return qstatxml

    def _qstatinfo(self, qstatxml):
        """
        Parameters
        ----------
        qstatxml : string
            The xml string of the 'qstat -xml' call.
        Returns
        -------
        queue_info : list
            A list of jobs in 'queue_info'. Jobs are dictionaries with both
            string keys and string names.
        job_info : list
            A list of jobs in 'job_info'.
        """
        x = xmltodict.parse(qstatxml)
        queue_info = []
        if x['job_info']['queue_info'] is not None:
            for job in x['job_info']['queue_info']['job_list']:
                queue_info.append(dict(job))
        job_info = []
        if x['job_info']['job_info'] is not None:
            for job in x['job_info']['job_info']['job_list']:
                job_info.append(dict(job))
        return queue_info, job_info
