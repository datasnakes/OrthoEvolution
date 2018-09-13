"""Other utilities optimized for this package/project."""
import contextlib
import pkg_resources
from importlib import import_module
import os
import psutil
from threading import Timer
from subprocess import run, CalledProcessError, PIPE
import pandas as pd
from pathlib import Path


class OtherUtils(object):
    def __init__(self):
        pass

    def splitlist(self, listname, basefilename, n):
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

    def formatlist(self, input_list):
        """Remove spaces from list items and turn those spaces into underscores."""
        output_list = []
        for item in input_list:
            item = str(item)
            item = item.replace(" ", "_")
            output_list.append(item)
            return output_list

    def makedirectory(self, path):
        """Creates path/parents and is compatible for python 3.4 and upwards."""
        exist_ok = True
        if not exist_ok and os.path.isdir(path):
            with contextlib.suppress(OSError):
                Path.mkdir(path, parents=True)

    # def set_paths(self, parent, **children):
    #     raise NotImplementedError("This function is being developed.")

    def csvtolist(self, csvfile, column_header='Organism'):
        """Turn column from csv file into a list."""
        file = pd.read_csv(csvfile)
        # Create a list name/variable and use list()
        listfromcolumn = list(file[column_header])

        return listfromcolumn

    def runcmd(self, command_string):
        """Run a command string.

        :param command string:
        """
        try:
            cmd = [command_string]  # this is the command
            # Shell MUST be True
            cmd_status = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, check=True)
        except CalledProcessError as err:
            errmsg = err.stderr.decode('utf-8')
            return errmsg
        else:
            if cmd_status.returncode == 0:  # Command was successful.
                # The cmd_status has stdout that must be decoded.
                cmd_stdout = cmd_status.stdout.decode('utf-8')
                return cmd_stdout
            else:  # Unsuccessful. Stdout will be '1'
                failmsg = '%s failed.' % command_string
                return failmsg

    # Determine if there are processes with this file opened (Linux)
    def has_handle(self, fpath):
        for proc in psutil.process_iter():
            try:
                for item in proc.open_files():
                    if fpath == item.path:
                        return True
            except Exception:
                pass
        return False

    # Safely open a file thats being accessed by multiple processes
    def safe_open(self, fpath, mode, iterations, cnt=0):
        if cnt == iterations:
            raise IOError
        unsafe = self.has_handle(fpath)
        if unsafe is True:
            cnt = cnt + 1
            self.safe_open(fpath=fpath, mode=mode, iterations=iterations, cnt=cnt)
        else:
            return open(file=fpath, mode=mode)


class PackageVersion(object):
    """Get the version of an installed python package."""
    def __init__(self, packagename):
        self.packagename = packagename
        self._getversion()

    def _getversion(self):
        import_module(self.packagename)
        version = pkg_resources.get_distribution(self.packagename).version
        print('Version %s of %s is installed.' % (version, self.packagename))


class FunctionRepeater(object):
    """Repeats a function every interval. Ref: https://tinyurl.com/yckgv8m2"""
    def __init__(self, interval, function, *args, **kwargs):
        self._timer = None
        self.function = function
        self.interval = interval
        self.args = args
        self.kwargs = kwargs
        self.is_running = False
        self.start()

    def _run(self):
        self.is_running = False
        self.start()
        self.function(*self.args, **self.kwargs)

    def start(self):
        if not self.is_running:
            self._timer = Timer(self.interval, self._run)
            self._timer.start()
            self.is_running = True

    def stop(self):
        self._timer.cancel()
        self.is_running = False
