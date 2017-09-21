"""Use python's multiprocessing module to create multiple processes and speed
up the completion of functions/classes."""

from multiprocessing import Pool, cpu_count
from time import time
import sys
from logzero import logger as log


class Multiprocess(object):
    """Use multiple processes with functions."""
    def __init__(self, num_procs, function, listinput):
        self.num_procs = num_procs  # Number of processes to run

        # Force user to use optimal number of processes
        if int(self.num_procs) > (int(cpu_count) - 1):
            sys.exit(ValueError)

        self.function = function
        self.listinput = listinput

    @classmethod
    def map2function(cls, num_procs, function, listinput):
        """Start a pool to run your function with a list."""
        time_secs = time()
        with Pool(processes=num_procs) as pool:
            pool.map(function, listinput)
            minutes = (time() - time_secs) / 60
            log.info("Took %s minutes to complete.", minutes)
