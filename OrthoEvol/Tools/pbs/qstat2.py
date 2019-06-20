import os
import csv
import yaml
import subprocess as sp
import pandas as pd
from dateutil import parser
from datetime import datetime
from pkg_resources import resource_filename
from collections import OrderedDict
from pathlib import Path
from OrthoEvol.utilities import FullUtilities
from OrthoEvol.Manager.config import yml
from OrthoEvol.Tools.logit import LogIt


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

        self.qstat_utils = FullUtilities()
        self.qstat_log = LogIt().default(logname="PBS - QSTAT", logfile=None)
        self._yaml_config = resource_filename(yml.__name__, 'qstat_dict.yml')
        self.cmd = cmd
        self.target_job = job
        self.outfile = outfile
        if not home:
            self.home = Path(os.getcwd()) / str(job).replace(".", "")
        else:
            self.home = Path(home)
        if not self.home.exists():
            self.home.mkdir(parents=True)

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
            # For this to work the infile must be an absolute path.
            if Path(infile).exists():
                self.configure_data_file(extra_data=infile)
            else:
                raise FileExistsError("The infile must be an absolute path.")

        # QSTAT data objects
        self.qstat_data = None
        self.qstat_dict = None
        self.static_data = None
        self.qstat_dataframe = None
        self.target_job_dict = None
        self.target_job_dataframe = None

    def configure_data_file(self, extra_data):
        """
        Configure the primary data file by appending extra data to it.  The csv header
        is only appended when the primary data file does not exist.

        :param extra_data:  The abolute path to a csv file that contains qstat data.
        :type extra_data:  str.
        """
        # Open infile containing data
        with open(extra_data, 'r') as _if:
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

    def get_qstat_output(self):
        """
        A function that calls qstat in a subprocess and stores the output
        in a class variable (qstat_data).  The data is the list returned from
        readlines().
        """
        try:
            proc = self.qstat_utils.system_cmd(self.cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True, encoding='utf-8',
                                               universal_newlines=False)
        except sp.CalledProcessError as err:
            self.qstat_log.error(err.stderr.decode('utf-8'))
        else:
            if proc.returncode == 0:
                self.qstat_data = proc.stdout.readlines()

    def qstat_to_dict(self, qstat_data):
        """
        The qstat parser takes the qstat data from the 'qstat -f' command and parses it
        into a ditionary.  It uses the qstat keywords found in the qstat yaml file.
        """
        mast_dict = OrderedDict()
        job_count = 0
        phrase_continuation_flag = None
        with open(self._yaml_config, 'r') as yf:
            qstat_keywords = yaml.load(yf)
        if isinstance(qstat_data, list):
            qstat_sentence = None
            continuation_phrase = ""
            qstat_phrase = ""
            prev_item = None
            for item in qstat_data:
                # If a new job is identified then create the nested dictionary
                if "Job Id" in item:
                    job_count += 1
                    _ = item.split(": ")
                    job_id_key = "%s" % _[1].replace("\r\n", "")
                    job_id_key = job_id_key.replace("\n", "")
                    mast_dict[job_id_key] = OrderedDict()
                    mast_dict[job_id_key]["Job_Id"] = job_id_key
                # The current line information is used to determine single or multi-lined parsing for the
                # previous line.
                # If the a new keyword is recognized, then parse the line.
                elif "    " in item and any(kw in item for kw in list(qstat_keywords["Job Id"].keys()) + list(qstat_keywords["Job Id"]["Variable_List"].keys())) or item == "\n":
                    item = item.replace("    ", "")
                    # Join the multi-lined "phrases" into one "sentence"
                    if phrase_continuation_flag is True:
                        qstat_sentence = qstat_phrase + continuation_phrase
                        phrase_continuation_flag = False
                    #  If the phrase is a single line
                    elif qstat_phrase == prev_item:
                        qstat_sentence = qstat_phrase
                    # Updates the qstat phrase unless it's in between lines or the end of file
                    if item != "\n":
                        qstat_phrase = item
                    else:
                        pass
                # If there is no keyword and tabbed whitespace is recognized, then the current line
                # is a continuation phrase for the most recent qstat keyword
                elif "\t" in item:
                    phrase_continuation_flag = True
                    continuation_phrase = continuation_phrase + item

                # For multi-line phrases format the qstat sentence
                if phrase_continuation_flag is False:
                    qstat_sentence = qstat_sentence.replace("\n\t", "").replace('\n', '')
                    continuation_phrase = ""
                    phrase_continuation_flag = None
                    qstat_list = qstat_sentence.split(" = ")
                    qstat_keyword = qstat_list[0]
                    qstat_value = qstat_list[1]
                    # The variable list is unique in that it can be split into a dictionary of
                    # environment variables
                    if qstat_keyword == "Variable_List":
                        qstat_value = qstat_value.split(",")
                        temp_dict = OrderedDict(var.split("=") for var in qstat_value)
                        for vl_key, vl_value in temp_dict.items():
                            if vl_value[0] == "/" or vl_value[0] == "\\":
                                vl_value = Path(vl_value)
                            temp_dict[vl_key] = vl_value
                        mast_dict[job_id_key]["Variable_List"] = temp_dict
                    # All of the other qstat keywords/sentences are basic key/value pairs
                    else:
                        mast_dict[job_id_key][qstat_keyword] = qstat_value

                # For single line Phrases
                elif qstat_sentence:
                    qstat_sentence = qstat_sentence.replace("\n\t", "")
                    qstat_list = qstat_sentence.split(" = ")
                    qstat_keyword = qstat_list[0]
                    qstat_value = qstat_list[1]
                    mast_dict[job_id_key][qstat_keyword] = qstat_value.replace('\n', '')
                qstat_sentence = None
                prev_item = item
        return mast_dict

    def filter_qstat_jobs(self, qstat_dict, target_job):
        """
        Filter out all of the qstat jobs other than the target job.  This sets up the
        class variables for the target dataframe and target dictionary.
        """
        target_job_dict = None
        for j in qstat_dict.keys():
            if qstat_dict[j]["Job_Id"] == target_job:
                target_job_dict = qstat_dict[j]
                break
        if not target_job_dict:
            raise ValueError("The target job does not exist in the qstat data provided.")
        else:
            return target_job_dict

    def filter_qstat_keywords(self, qstat_dict, static_flag=False, python_datetime=datetime.now()):
        """
        This function takes a qstat dictionary and returns a dictionary that contains "dynamic" data that can be plotted
        or "static" data related to the jobs.
        :param qstat_dict:
        :type qstat_dict:
        :param static_flag:
        :type static_flag:
        :param python_datetime:
        :type python_datetime:
        :return:
        :rtype:
        """
        data_dict = OrderedDict()
        if static_flag:
            for job in qstat_dict.keys():
                data_dict[job] = OrderedDict()
                for keyword in qstat_dict[job].keys():
                    if keyword in self.__static_kw:
                        if keyword in self.__job_time_kw:
                            data_dict[job][keyword] = str(parser.parse(qstat_dict[job][keyword]))
                        else:
                            data_dict[job][keyword] = qstat_dict[job][keyword]
            _data = data_dict
        else:
            for job in qstat_dict.keys():
                data_dict[job] = OrderedDict()
                # Store the python datetime
                data_dict[job]["datetime"] = [python_datetime]
                for keyword in qstat_dict[job].keys():
                    # Store all of the dynamic data so that it can be converted to a dataframe.
                    if keyword in self.__dynamic_kw:
                        data_dict[job][keyword] = [qstat_dict[job][keyword]]

            # Finish converting dynamic data into a pandas dataframe
            if len(qstat_dict.keys()) > 0:
                data_frame = pd.DataFrame.from_dict(data_dict)
            else:
                for job in qstat_dict.keys():
                    data_frame = pd.DataFrame.from_dict(dict(qstat_dict[job]))
            _data = data_frame
        return _data
