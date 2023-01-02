"""Main logging class to make logging easier."""
import os
import sys
from logzero import setup_logger, LogFormatter, logging, colors
from logging import CRITICAL, ERROR, WARNING, INFO, DEBUG

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

        # Add custom colors for CRITICAL and DEBUG
        self._COLORS = {DEBUG: colors.Fore.LIGHTBLUE_EX,
        INFO: colors.Fore.GREEN,
        WARNING: colors.Fore.YELLOW,
        ERROR: colors.Fore.RED,
        CRITICAL: colors.Fore.LIGHTRED_EX
        }
        
        self._formatter = LogFormatter(fmt=self._log_format,
                                       datefmt=self._date_format,
                                       colors=self._COLORS)

        self.logging = logging

    def default(self, logname, logfile):
        """Create a log handler using default formatting.

        :param logname: Name of the log.
        :param logfile: Path to the logfile. Can be `None`.
        :return: A default_log object is returned.
        """

        default_log = setup_logger(name=logname.upper(), logfile=logfile,
                                   level=logging.DEBUG,
                                   formatter=self._formatter)

        assert isinstance(default_log, object)
        return default_log

    @classmethod
    def custom(cls, logname, logfile, level, fmt='default'):
        """Create a log handler or logger.

        :param logname: Name of the log.
        :param logfile: Path to the logfile. Can be `None`
        :param level: Logging level.
        :param fmt: Logging format (fmt is default).
        :return: A custom_log object is returned.
        """

        if fmt == 'default':
            fmt = '[%(levelname)-2s - %(name)s]: %(message)s'
        elif fmt == 'custom':
            # TODO Allow customization
            raise NotImplementedError('Feature not integrated yet!')
        elif fmt != 'default' or 'custom':
            raise Exception('User did not provide a format for the logger.')

        custom_log = setup_logger(name=logname, logfile=logfile,
                                  level=level, formatter=fmt)
        return custom_log

    def deletelog(self, logfile):
        """Delete the log file.

        :param logfile: Name of the logfile to be deleted.
        """

        self.shutdown()
        if os.path.isfile(logfile):
            os.remove(logfile)

    def shutdown(self):
        """Shutdown the log handlers.

        Windows won't just delete a log (linux distros will).
        It must be shutdown.
        """

        # TIP https://www.programcreek.com/python/example/3517/logging.shutdown
        self.logging.shutdown()
