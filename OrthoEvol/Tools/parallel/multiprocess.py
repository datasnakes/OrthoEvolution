"""Use python's multiprocessing module to create multiple processes and speed
up the completion of functions/classes."""

from multiprocessing import Pool, cpu_count, get_logger
from time import time
import logging

import logzero


class Multiprocess(object):
    """Use multiple processes with a function."""

    cpus = cpu_count()
    num_procs = cpus - 1

    def __init__(self):
        pass

    @staticmethod
    def _logger():
        """Add the multiprocessing module's logger.

        :return: Returns a multiprocessing logger.
        """

        multiprocess_handler = get_logger()
        multiprocess_handler = logging.StreamHandler()
        multiprocess_handler.setLevel(logging.ERROR)
        multiprocess_handler.setFormatter(logzero.LogFormatter(color=True))

        # Attach it to the logzero default logger
        logzero.logger.addHandler(multiprocess_handler)
        logger = logzero.logger
        return logger

    def map2function(self, function, iterable):
        """Start a pool to run your function with a list.

        :param function: Input a python function.
        :param iterable: Input a list or dictionary to map to the function.
        """

        log = self._logger()  # Start the logger
        time_secs = time()
        with Pool(processes=self.num_procs) as pool:
            pool.map(function, iterable)
            minutes = (time() - time_secs) / 60
        log.info("Took %s minutes to complete.", minutes)
        logging.shutdown()  # Shutdown the logger.
