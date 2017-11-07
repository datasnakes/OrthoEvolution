"""Collection of tools for using PBS, a job scheduler for high-performance
computing environments. The command is usually `qsub <options>` on most
systems.
"""
import random
import string
from string import Template


def basejobids(length, name='submit'):
    """"Create base job id and name."""
    base_id = randomid(length=length)
    base = name + "_{0}".format(base_id)

    return base_id, base


def import_temp(filepath):
    """Import the script or file that you need a template of and that has
    temp strings.
    """
    file_temp = open(filepath, 'r')
    file_str = file_temp.read()
    file_temp.close()

    file_temp = Template(file_str)
    return file_temp


def file2str(filepath):
    """Turn the contents of a file (python file) into a string."""
    file_temp = open(filepath, 'r')
    file_str = file_temp.read()
    return file_str


def randomid(length=5):
    """Generate a random ID of 5 characters to append to qsub job name."""
    return ''.join(random.sample(string.ascii_letters + string.digits, length))


def writecodefile(filename, code, language):
    """Create a python file and write the code to it."""
    if language == 'python':
        with open(filename + '.py', 'w') as pyfile:
            pyfile.write(code)
            pyfile.close()

    elif language == 'bash':
        with open(filename + '.sh', 'w') as bashfile:
            bashfile.write(code)
            bashfile.close()

    elif language == 'R' or 'r':
        with open(filename + '.R', 'w') as rfile:
            rfile.write(code)
            rfile.close()
    else:
        raise NotImplementedError('%s is unsupported.' % language)
