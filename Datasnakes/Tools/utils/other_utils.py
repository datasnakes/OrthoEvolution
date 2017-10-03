"""Other utilities optimized for this package/project."""
import pandas as pd
import contextlib
from pathlib import Path
import pkg_resources
from importlib import import_module
import os


def splitlist(listname, basefilename, n):
    """Split a long list into chunks and save chunks as a text file."""
    # Split the list into chunks
    chunks = [listname[x:x + n] for x in range(0, len(listname), n)]
    list_group = []
    num_lists = len(chunks)

    # Name and save the lists
    for chunk, num in zip(chunks, range(0, num_lists)):
        listdf = pd.DataFrame(chunk)
        n = basefilename + '_list_' + str(num)
        listdf.to_csv(n + ".txt", index=False, header=None)
        list_group.append(n)
    return list_group


def formatlist(input_list):
    """Remove spaces from list items and turn those spaces into underscores."""
    output_list = []
    for item in input_list:
        item = str(item)
        item = item.replace(" ", "_")
        output_list.append(item)
        return output_list


def makedirectory(path):
    """Creates path/parents and is compatible for python 3.4 and upwards."""
    exist_ok = True
    if not exist_ok and os.path.isdir(path):
        with contextlib.suppress(OSError):
            Path.mkdir(path, parents=True)


class PackageVersion(object):
    """Get the version of an installed python package."""
    def __init__(self, packagename):
        self.packagename = packagename
        self._getversion()

    def _getversion(self):
        import_module(self.packagename)
        version = pkg_resources.get_distribution(self.packagename).version
        print('Version %s of %s is installed.' % (version, self.packagename))

def set_paths(parent, **children):
    raise NotImplementedError("This function is being developed.")