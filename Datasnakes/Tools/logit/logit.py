"""Main logging class to make logging easier."""
import os
import sys
from logzero import setup_logger, LogFormatter, logging


class LogIt(object):
    """LogIt makes logging easier by creating easy loggers."""
    def __init__(self):
        """Initialize the logger format based on system platform."""
        # Set the different formats based on user's platform
        if sys.platform == 'win32':
            self._archive_format = '%m-%d-%Y_%I-%M-%p'
        elif sys.platform == 'linux':
            self._archive_format = '%m-%d-%Y@%I:%M:%S-%p'

        self._date_format = '%b-%d-%Y at %I:%M:%S %p'  # Used to add date
        self._log_format = ("%(color)s[%(levelname)s | %(name)s] [%(asctime)s | "
                            "%(module)s - line %(lineno)d]:%(end_color)s %(message)s")
        self._formatter = LogFormatter(fmt=self._log_format,
                                       datefmt=self._date_format)
        self.logging = logging

    def default(self, logname, logfile):
        """Create a log handler using default formatting."""
        default_log = setup_logger(name=logname.upper(), logfile=logfile,
                                   level=logging.DEBUG, formatter=self._formatter)
        return default_log

    @classmethod
    def custom(cls, logname, logfile, level, fmt='default'):
        """Create a log handler or logger."""
        if fmt is 'default':
            fmt = '[%(levelname)-2s - %(name)s]: %(message)s'
        elif fmt is 'custom':
            # TODO Allow customization
            raise NotImplementedError('Feature not integrated yet!')
        elif fmt is not 'default' or 'custom':
            raise Exception('User did not provide a format for the logger.')

        custom_log = setup_logger(name=logname, logfile=logfile,
                                  level=level, formatter=fmt)
        return custom_log

    def deletelog(self, logfile):
        """Delete the log file."""
        self.shutdown()
        # TODO Use contextlib here; See makedirectory function
        if os.path.isfile(logfile):
            os.remove(logfile)

    def shutdown(self):
        """Shutdown the log handlers."""
        # HINT https://www.programcreek.com/python/example/3517/logging.shutdown
        self.logging.shutdown()
        # Windows won't just delete a log (linux distros will).
        # It must be shutdown.
