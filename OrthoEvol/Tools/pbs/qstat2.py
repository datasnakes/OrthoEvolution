

class BaseQstat(object):
    """
    The BaseQstat class processes the output from the pbs command 'qstat'.  It
    specifically parses output from 'qstat -f', which displays a full status report
    for all of the jobs in the queue.  Because each line of output per job consists of
    attribute_names and values for those attributes, the qstat data is parsed into a
    dictionary.  The qstat data is then converted to csv format and saved in a .csv file.
    """
    def __init__(self):
        pass