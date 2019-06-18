from pathlib import Path
import os
import csv


class BaseQstat(object):
    # Static qstat Keywords

    __misc_kw = ["Checkpoint", "Error_Path", "exec_host", "exec_vnode", "Hold_Types", "Join_Path",
                             "Keep_Files", "Mail_Points", "Output_Path", "Rerunable", "Resource_List.mpiprocs",
                             "Resource_List.ncpus", "Resource_List.nodect", "Resource_List.nodes",
                             "Resource_List.place", "Resource_List.select", "jobdir", "Variable_List", "umask",
                             "project", "Submit_arguments"]
    __job_limits_kw = ["ctime", "etime", "qtime", "stime", "mtime", "Resource_List.walltime", "Resource_List.cput",
                           "Resource_List.mem"]
    __job_time_kw = ["ctime", "etime", "qtime", "stime", "mtime"]
    __job_info_kw = ["Job_Id", "Job_Name", "Job_Owner", "queue", "server", "session_id"]
    __static_kw = __job_info_kw + __job_limits_kw + __misc_kw
    # Dynamic qstat Keywords
    __misc_data_kw = ["job_state", "Priority", "substate", "comment", "run_count"]
    __job_data_kw = ["resources_used.cpupercent", "resources_used.cput", "resources_used.mem",
                            "resources_used.vmem", "resources_used.walltime", "resources_used.ncpus"]
    __dynamic_kw = __job_data_kw + __misc_data_kw
    # All Keywords
    __keywords = __static_kw + __dynamic_kw
    # Metadata options
    __metadata_dict = {"environment variables": "Variable_List",
                       "plot": {"limits": __job_limits_kw,
                                "info": __job_info_kw,
                                "data": __job_data_kw},
                       "all": __keywords
                       }

    def __init__(self, job, infile=None, outfile=None, home=None, cmd="qstat -f"):
        """
        The BaseQstat class processes the output from the pbs command 'qstat'.  It
        specifically parses output from 'qstat -f', which displays a full status report
        for all of the jobs in the queue.  Because each line of output per job consists of
        attribute_names and values for those attributes, the qstat data is parsed into a
        dictionary.  The qstat data is then converted to csv format and saved in a .csv file.
        The Base class only process one job, and it only gathers data on one point in time.

        :param job:  The name of the job to analyze.
        :type job:  str.
        :param infile:  The input file and the output file are used in tandem to determine the
        data file that will be used.  If only one of these values are given (infile/outfile), then
        it will be used as the data file.  If neither of these values are given, then a default file
        ("job_data.csv") will be used.  If both are given, then the infile data is appended to the outfile,
        which is used as the data file.
        :type infile:  str.
        :param outfile:  See infile.
        :type outfile: str.
        :param home:  An absolute path to the director where qstat data will be stored.
        :type home:  str.
        :param cmd:  The qstat command used to produce the job status report.
        :type cmd:  str.
        """

        self.job = job
        self.outfile = outfile
        if not home:
            self.home = Path(os.getcwd()) / str(job).replace(".", "")
        else:
            self.home = home

        # Use infile as data file whether it exists or not
        if infile is not None and outfile is None:
            self.data_file = self.home / infile
        # Use outfile as data file whether it exists or not
        elif infile is None and outfile is not None:
            self.data_file = self.home / outfile
        # Use a default file name like "job_data.csv"
        elif infile is None and outfile is None:
            self.data_file = self.home / "job_data.csv"
        # Use the outfile as the data file whether it exists or not,
        # but also append the infile data to it.
        elif infile is not None and outfile is not None:
            self.data_file = self.home / outfile
            self.configure_data_file(infile=infile)

        self.cmd = cmd

    def configure_data_file(self, infile):
        # For this to work the infile must be an absolute path.
        if Path(infile).exists():
            # Open infile containing data
            with open(infile, 'r') as _if:
                in_data = csv.reader(_if, delimiter=",")
                line_count = 0
                header_flag = Path(self.data_file).exists()
                with open(self.data_file, 'a') as _df:
                    out_data = csv.writer(_df, delimiter=",")
                    for row in in_data:
                        if line_count == 0:
                            # Write infile header if data file doesn't exist
                            if not header_flag:
                                out_data.writerow(row)
                            line_count += 1
                        else:
                            out_data.writerow(row)
