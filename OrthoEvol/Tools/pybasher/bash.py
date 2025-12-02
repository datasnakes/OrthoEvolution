"""A user-friendly Bash module for Python.

Built upon Alex Couper's `bash` package. (https://github.com/alexcouper/bash)
"""
# TODO-SDH It may be helpful to use psutil for pybasher.
# TODO-SDH Finish pybasher ASAP
# TODO-SDH Look at some examples for pybasher.
# TODO-SDH Update the README for pybasher
import platform
import sys

from OrthoEvol.utilities import FullUtilities
# Set up run cmd function
runcmd = FullUtilities().system_cmd

SUBPROCESS_HAS_TIMEOUT = True

if "windows" in platform.system().lower():
    raise ImportError("PyBasher is currently only supported on linux and osx.")
else:
    from subprocess import PIPE, Popen

    if sys.version_info < (3, 0):
        try:
            from subprocess32 import PIPE, Popen
        except ImportError:
            # You haven't got subprocess32 installed. If you're running 2.X this
            # will mean you don't have access to things like timeout
            SUBPROCESS_HAS_TIMEOUT = False

#import os
#import configparser
# TODO-SDH use a config file to load/use a list or group of common commands.


class BaseBash(object):
    """Utilize bash commands within python."""

    def __init__(self):
        """Initialize the call as well as standard error and output."""
        # TODO-SDH Test if this is working.
        self.process = None
        self.stdout = None

    def _bash(self, cmd, env=None, stdout=PIPE, stderr=PIPE, timeout=None, _sync=True):
        """Use subprocess to run bash commands.

        :param cmd: The bash command to be run.
        :param env:  (Default value = None)
        :param stdout:  (Default value = PIPE)
        :param stderr:  (Default value = PIPE)
        :param timeout:  (Default value = None)
        :param _sync:  (Default value = True)
        """

        self.process = Popen(cmd, shell=True, stdout=stdout, stdin=PIPE,
                             stderr=stderr, env=env)
        if _sync:
            self._sync(timeout)
        return self

    def _sync(self, timeout=None):
        """Ensure function is run.

        :param timeout:  (Default value = None)
        """

        kwargs = {'input': self.stdout}
        if timeout:
            kwargs['timeout'] = timeout
            if not SUBPROCESS_HAS_TIMEOUT:
                raise ValueError(
                    "Timeout given but subprocess doesn't support it. "
                    "Install subprocess32 and try again."
                )
        self.stdout, self.stderr = self.process.communicate(**kwargs)
        self.code = self.process.returncode
        return self

    def __repr__(self):
        return self._value()

    def __unicode__(self):
        return self._value()

    def __str__(self):
        return self._value()

    def __nonzero__(self):
        return self.__bool__()

    def __bool__(self):
        return bool(self.value())

    def _value(self):
        if self.stdout:
            return self.stdout.strip().decode(encoding='UTF-8')
        return ''


class PyBasher(BaseBash):
    """Common bash commands."""

    def __init__(self):
        super().__init__()

    def cp(self, source, destination):
        """Copy file."""

        cmd = f"cp {source} {destination}"
        self._bash(cmd)

    def mv(self, source, destination):
        """Move file."""

        cmd = f"mv {source} {destination}"
        self._bash(cmd)

    def rm(self, path):
        """Delete file."""

        cmd = f"rm {path}"
        self._bash(cmd)

    def mkdir(self, path):
        """Create directory."""

        cmd = f"mkdir {path}"
        self._bash(cmd)

    def rmdir(self, path):
        """Delete empty directory."""

        cmd = f"rmdir {path}"
        self._bash(cmd)

    def touch(self, path):
        """Create empty file."""

        cmd = f"touch {path}"
        self._bash(cmd)

    def ls(self, path):
        """List directory contents."""

        cmd = f"ls {path}"
        self._bash(cmd)

    def cat(self, path):
        """Display contents of file."""

        cmd = f"cat {path}"
        self._bash(cmd)

    def pwd(self):
        """Print current working directory."""

        cmd = "pwd"
        self._bash(cmd)

    def cd(self, path):
        """Change current working directory."""

        cmd = f"cd {path}"
        self._bash(cmd)

    def grep(self, pattern, path):
        """Search for pattern in file."""

        cmd = f"grep {pattern} {path}"
        self._bash(cmd)

    def chmod(self, permissions, path):
        """Change file permissions."""

        cmd = f"chmod {permissions} {path}"
        self._bash(cmd)

    def chown(self, owner, path):
        """Change file owner."""

        cmd = f"chown {owner} {path}"
        self._bash(cmd)

    def find(self, path, pattern):
        """Search for files matching pattern in path."""

        cmd = f"find {path} -name {pattern}"
        self._bash(cmd)

    def tar(self, action, options, archive, files):
        """Create or extract tar archive."""

        cmd = f"tar {action} {options} {archive} {files}"
        self._bash(cmd)

    def unzip(self, archive, destination):
        """Extract zip archive."""

        cmd = f"unzip {archive} -d {destination}"
        self._bash(cmd)

    def zip(self, options, archive, files):
        """Create zip archive."""

        cmd = f"zip {options} {archive} {files}"
        self._bash(cmd)
        
    def du(self, path):
        """Display disk usage statistics."""

        cmd = f"du {path}"
        self._bash(cmd)

    def df(self, path):
        """Display free disk space."""

        cmd = f"df {path}"
        self._bash(cmd)

    def top(self, options):
        """Display top-running processes."""

        cmd = f"top {options}"
        self._bash(cmd)

    def ps(self, options):
        """Display running processes."""

        cmd = f"ps {options}"
        self._bash(cmd)

    def kill(self, signal, process_id):
        """Send signal to process."""

        cmd = f"kill {signal} {process_id}"
        self._bash(cmd)

    def man(self, command):
        """Display manual for command."""

        cmd = f"man {command}"
        self._bash(cmd)

    def info(self, command):
        """Display information for command."""

        cmd = f"info {command}"
        self._bash(cmd)

    def history(self, options):
        """Display command history."""

        cmd = f"history {options}"
        self._bash(cmd)