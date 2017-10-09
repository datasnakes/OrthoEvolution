"""Use python's multiprocessing module to create multiple processes and speed
up the completion of functions/classes."""

from multiprocessing import Pool, cpu_count, get_logger
from time import time
import logging

import logzero


cpus = cpu_count()
num_procs = cpus - 1


def _logger():
    """Add the multiprocessing module's logger."""
    multiprocess_handler = get_logger()
    multiprocess_handler = logging.StreamHandler()
    multiprocess_handler.setLevel(logging.ERROR)
    multiprocess_handler.setFormatter(logzero.LogFormatter(color=True))

    # Attach it to the logzero default logger
    logzero.logger.addHandler(multiprocess_handler)
    logger = logzero.logger
    return logger


def map2function(function, listinput):
    """Start a pool to run your function with a list."""
    log = _logger()  # Start the logger
    time_secs = time()
    with Pool(processes=num_procs) as pool:
        pool.map(function, listinput)
        minutes = (time() - time_secs) / 60
    log.info("Took %s minutes to complete.", minutes)
    logging.shutdown()  # Shutdown the logger
