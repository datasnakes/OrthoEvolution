"""Use python's multiprocessing module to create multiple processes and speed
up the completion of functions/classes."""

from multiprocessing import Pool, cpu_count, get_logger
from time import time
import logging
import logzero


class Multiprocess(object):
    """Use multiple processes with a function."""
    def __init__(self, num_procs, function, listinput):
        self.num_procs = num_procs  # Number of processes to run

        # Force user to use optimal number of processes
        if int(self.num_procs) > (int(cpu_count) - 1):
            raise ValueError('Too many processes assigned.')

        self.function = function  # User function that takes 1 item or list
        self.listinput = listinput  # List to map to a function
        self.iterable = self.listinput
        self.log = self._logger()

    def _logger(self):
        """Add the multiprocessing module's logger."""
        multiprocess_handler = get_logger()
        multiprocess_handler = logging.StreamHandler()
        multiprocess_handler.setLevel(logging.DEBUG)
        multiprocess_handler.setFormatter(logzero.LogFormatter(color=True))

        # Attach it to the logzero default logger
        logzero.logger.addHandler(multiprocess_handler)
        logger = logzero.logger
        return logger

    def map2function(self):
        """Start a pool to run your function with a list."""
        time_secs = time()
        with Pool(processes=self.num_procs) as pool:
            pool.map(self.function, self.iterable)
            minutes = (time() - time_secs) / 60
            self.log.info("Took %s minutes to complete.", minutes)
