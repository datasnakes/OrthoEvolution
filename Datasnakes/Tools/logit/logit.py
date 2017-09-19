"""Main logging class to make logging easier."""
from logzero import setup_logger, LogFormatter, logging
import os
import sys


class LogIt(object):
    """LogIt makes logging easier by creating easy loggers."""
    def __init__(self):
        """Initialize the logger format based on system platform."""
        # Set the different formats based on user's platform
        if sys.platform == 'win32':
            self.archive_format = '%m-%d-%Y_%I-%M-%p'
        elif sys.platform == 'linux':
            self.archive_format = '%m-%d-%Y@%I:%M:%S-%p'

        self.date_format = '%b-%d-%Y at %I:%M:%S %p'  # Used to add date
        self.log_format = '%(color)s[%(levelname)s | %(name)s] [%(asctime)s | %(module)s - line %(lineno)d]:%(end_color)s %(message)s'
        self.formatter = LogFormatter(fmt=self.log_format,
                                      datefmt=self.date_format)

        # TODO Allow customization
#        if not default:
#            self.custom_log = self.custom(logname, logfile, level, fmt)

    def default(self, logname, logfile):
        """Create a log handler using default formatting."""
        default_log = setup_logger(name=logname.upper(), logfile=logfile,
                                   level=logging.DEBUG, formatter=self.formatter)
        return default_log

    def custom(self, logname, logfile, level, fmt='default'):
        """Create a log handler or logger."""
        if fmt is 'default':
            fmt = '[%(levelname)-2s - %(name)s]: %(message)s'
        elif fmt is 'custom':
            fmt = self.formatcreator()
        elif fmt is not 'default' or 'custom':
            raise Exception('User did not provide a format for the logger.')

        custom_log = setup_logger(name=logname, logfile=logfile, level=level,
                                  formatter=fmt)
        return custom_log

    def deletelog(self, logfile):
        """Delete the log file."""
        self.shutdown()
        # TODO Use contextlib here
        if os.path.exists(logfile):
            os.remove(logfile)

    def shutdown(self):
        """Shutdown the log handlers."""
        # HINT https://www.programcreek.com/python/example/3517/logging.shutdown
        logging.shutdown()
        # Windows won't just delete a log (linux distros will).
        # It must be shutdown.