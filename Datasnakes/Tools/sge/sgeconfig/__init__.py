"""Configuration dictionaries for SGEJob class."""
import getpass
from datetime import datetime as d
import sys
from Datasnakes.Tools.sge import randomid

if sys.version_info.major < 3:
    raise NotImplementedError('This is not designed for the python version in your \
                              path')
elif sys.version_info.major >= 3 and sys.version_info.minor < 5:
    raise NotImplementedError('This is not designed for the python version in your \
                              path')

_format1 = '%a %b %d %I:%M:%S %p %Y'

_jobname = 'orthoevol_{}'.format(randomid(length=4))

__DEFAULT__ = {
            'author': getpass.getuser().upper(),
            'description': 'This is a default pbs job.',
            'date': d.now().strftime(_format1),
            'PBS_JOBID': '${PBS_JOBID}',
            'PBS_O_WORKDIR': '${PBS_O_WORKDIR}',
            'proj_name': 'Datasnakes-Orthologs',
            'select': 3,
            'memgb': 6,
            'cput': '{cput}',
            'wt': '{walltime}',
            'job_name': _jobname,
            'outfile': _jobname + '.o',
            'errfile': _jobname + '.o',
            'script': '{script_name}',
            'log_name': _jobname + '.log',
            'cmd': 'python3 {script_name}.py',
            'email': 'n/a'
             }


__CUSTOM__ = {
            'author': '{author}',
            'description': '{description}',
            'date': d.now().strftime(_format1),
            'PBS_JOBID': '${PBS_JOBID}',
            'PBS_O_WORKDIR': '${PBS_O_WORKDIR}',
            'proj_name': '{project}',
            'select': '{int(select)}',
            'memgb': '{gigabytes}',
            'cput': '{cput}',
            'wt': '{walltime}',
            'job_name': '{jobname}',
            'outfile': '{outfile}',
            'errfile': '{errfile}',
            'script': '{script_name}',
            'log_name': '{logname}',
            'cmd': 'python3 {script_name}.py',
            'email': '{email}'
            }
