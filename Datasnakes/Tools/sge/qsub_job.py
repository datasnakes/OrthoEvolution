from Datasnakes.Tools.sge import Qsubutils


class BaseJob(object):
    """Create a class for simple jobs."""
    def __init__(self):
        qsubutils = Qsubutils()
        self.qsubutils = qsubutils


class Qsubjob(object):
    """Create multiple jobs & scripts for each job to run based on
    splitting a list into chunks.
    """
    # TODO-SDH Create a more simplified process/algorithm.
