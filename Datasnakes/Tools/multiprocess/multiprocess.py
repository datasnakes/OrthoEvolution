"""Use python's multiprocessing module to create multiple process and speed up
the completion of classes."""
from logzero import logger as log
from multiprocessing import Pool, cpu_count
from time import time
import sys


class Multiprocess(object):
    def __init__(self, n, function, listinput):
        self.n = n  # Number of processes to run

        # Force user to use optimal number of processes
        if int(self.n) > (int(cpu_count) - 1):
            sys.exit(ValueError)

        self.function = function
        self.listinput = listinput

    def main(self):
        """Start a pool to run your function with a list."""
        ts = time()
        with Pool(processes=self.n) as p:
            p.map(self.function, self.listinput)
            minutes = (time() - ts) / 60
            log.info("Took {} minutes to complete".format(minutes))
