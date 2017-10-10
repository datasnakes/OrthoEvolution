from Datasnakes.Tools.sge import QsubUtils


class MultiJobber(QsubUtils):
    """Create multiple jobs & scripts for each job to run based on
    splitting a list into chunks.
    """
    # TODO-SDH This needs testing.
    # TODO-SDH Create a more simplified process/algorithm.
