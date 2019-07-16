import asyncio
import os
import csv
import yaml
import sys
import subprocess as sp
import pandas as pd
import plotly.graph_objs as go
import plotly
from dateutil import parser
from datetime import datetime
from time import sleep
from pkg_resources import resource_filename
from collections import OrderedDict
from pathlib import Path
from OrthoEvol.utilities import FullUtilities
from OrthoEvol.Manager.config import yml
from OrthoEvol.Tools.logit import LogIt


class TargetJobKeyError(KeyError):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


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
                self.configure_data_file(file=self.data_file, extra_data=infile)
            else:
                raise FileExistsError("The infile must be an absolute path.")
        self.info_file = self.data_file.parent / (str(self.data_file.stem) + '.yml')

        # QSTAT data objects
        self.qstat_data = None
        self.qstat_dict = None
        self.job_dict = None
        self.static_dict = None
        self.job_dataframe = None

    def configure_data_file(self, file, extra_data):
        """
        Configure the primary data file by appending extra data to it.  The csv header
        is only appended when the primary data file does not exist.

        :param file:  The absolute path to a csv file that will be updated.
        :type file:  str.
        :param extra_data:  The abolute path to a csv file that contains extra qstat data.
        :type extra_data:  str.
        """
        data_file = Path(file)
        # Open infile containing data
        with open(extra_data, 'r') as _if:
            in_data = csv.reader(_if, delimiter=",")
            line_count = 0
            header_flag = data_file.exists()
            with open(str(data_file), 'a') as _df:
                out_data = csv.writer(_df, delimiter=",")
                for row in in_data:
                    if line_count == 0:
                        # Write infile header if data file doesn't exist
                        if not header_flag:
                            out_data.writerow(row)
                        line_count += 1
                    else:
                        out_data.writerow(row)

    def run_qstat(self, csv_flag=True, sqlite_flag=False):
        """
        This method runs the qstat command, generates qstat data, parses it into various formats,
        and saves the data if desired.

        :param csv_flag:  A flag that determines if the data is saved in a csv file.
        :type csv_flag:  bool.
        :param sqlite_flag:  A flag that determines if the data is saved in a sqlite database.
        :type sqlite_flag:  bool.
        """
        # Get raw qstat data
        self.qstat_data = self.qstat_output(cmd=self.cmd)
        # Convert raw data to nested dictionary
        self.qstat_dict = self.to_dict(qstat_data=self.qstat_data)
        # Isolate data for target PBS job
        self.job_dict = self.target_data(qstat_dict=self.qstat_dict, target_job=self.target_job)
        # Isolate static data for target PBS job
        self.static_dict = self.static_data(qstat_dict=self.qstat_dict, target_job=self.target_job)
        # Create a pandas dataframe for target PBS job, formatted for creating a CSV file.
        self.job_dataframe = self.to_dataframe(qstat_dict=self.qstat_dict, target_job=self.target_job)
        if csv_flag:
            self.to_csv(file=self.data_file, qstat_dict=self.qstat_dict, target_job=self.target_job)
            self.static_data_to_yaml(file=self.info_file, qstat_dict=self.qstat_dict, target_job=self.target_job)
        if sqlite_flag:
            self.to_sqlite()

    def qstat_output(self, cmd):
        """
        A function that calls qstat via subprocess.  The data is the list returned from
        readlines().

        :param cmd:  The qstat command used to generate qstat data.  This is usually
        'qstat -f'
        :type cmd:  str.
        :return:  Output generated and read from the qstat command.
        :rtype:  list.
        """
        try:
            proc = self.qstat_utils.system_cmd(cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True,
                                               universal_newlines=False)
        except sp.CalledProcessError as err:
            self.qstat_log.error(err.stderr.decode('utf-8'))
        else:
            if proc.returncode == 0:
                qstat_data = proc.stdout.readlines()
                return qstat_data

    def to_dict(self, qstat_data):
        """
        The qstat parser takes the qstat data from the 'qstat -f' command and parses it
        into a ditionary.  It uses the qstat keywords found in the qstat yaml file.

        :param qstat_data:  Output data generated from the qstat command.
        :type qstat_data:  list.
        :return:  A nest dictionary that uses JobId's as keys.
        :rtype:  dict.
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

    def target_data(self, qstat_dict, target_job):
        """
        Filter out all of the qstat jobs other than the target job.  This sets up the
        class variables for the target dataframe and target dictionary.

                :param qstat_dict:  Qstat data that has been parsed into a dictionary.
        :type qstat_dict:  dict.
        :param target_job:  The target job that's being analyzed.
        :type target_job:  str.
        :return:  A dictionary that contains all of the target job's static and dynamic data.
        :rtype:  dict.
        """
        if target_job not in qstat_dict.keys():
            raise ValueError("The target job does not exist in the qstat data provided.")
        target_job_dict = qstat_dict[target_job]
        return target_job_dict

    def static_data(self, qstat_dict, target_job):
        """
        This function takes a qstat dictionary and returns a dictionary that contains
        "static" data related to the job of interest.

        :param qstat_dict:  Qstat data that has been parsed into a dictionary.
        :type qstat_dict:  dict.
        :param target_job:  The target job that's being analyzed.
        :type target_job:  str.
        :return:  A dictionary that contains the static data of the target job.
        :rtype:  dict.
        """
        if target_job not in qstat_dict.keys():
            raise ValueError("The target job does not exist in the qstat data provided.")

        data_dict = OrderedDict()
        data_dict[target_job] = OrderedDict()
        for keyword in qstat_dict[target_job].keys():
            if keyword in self.__static_kw:
                if keyword in self.__job_time_kw:
                    data_dict[target_job][keyword] = str(parser.parse(qstat_dict[target_job][keyword]))
                else:
                    data_dict[target_job][keyword] = qstat_dict[target_job][keyword]
        return data_dict

    def static_data_to_yaml(self, file, qstat_dict, target_job, overwrite=False):
        """
        This method saves the static data to a YAML file.

        :param file:  The absolute path of the yaml file.
        :type file:  str.
        :param qstat_dict:  Qstat data that has been parsed into a dictionary.
        :type qstat_dict:  dict.
        :param target_job:  The target job that's being analyzed.
        :type target_job:  str.
        :param overwrite:  A flag to determine weather the yaml file will be overwritten.
        :type overwrite:  bool.
        """
        data_file = Path(file)
        static_data = self.static_data(qstat_dict=qstat_dict, target_job=target_job)
        if data_file.is_file():
            if not overwrite:
                with open(data_file, 'a') as _f:
                    yaml.dump(static_data, _f, default_flow_style=False)
            else:
                with open(data_file, 'w') as _f:
                    yaml.dump(static_data, _f, default_flow_style=False)
        else:
            with open(data_file, 'w') as _f:
                yaml.dump(static_data, _f, default_flow_style=False)

    def to_dataframe(self, qstat_dict, target_job, python_datetime=datetime.now()):
        """
        Convert a qstat dictionary to a dataframe that contains data, which can be
        plotted.

        :param qstat_dict:  Qstat data that has been parsed into a dictionary.
        :type qstat_dict:  dict.
        :param target_job:  The target job that's being analyzed.
        :type target_job:  str.
        :param python_datetime:  The date and time that the qstat data was collected.
        :type python_datetime:  str.
        :return:  A dataframe that contains dynamic data.
        :rtype:  pd.DataFrame
        """
        if target_job not in qstat_dict.keys():
            raise ValueError("The target job does not exist in the qstat data provided.")

        data_dict = OrderedDict()
        # Store the python datetime
        data_dict["datetime"] = [python_datetime]
        data_dict["Job_Id"] = [target_job]
        for keyword in qstat_dict[target_job].keys():
            # Store all of the dynamic data so that it can be converted to a dataframe.
            if keyword in self.__dynamic_kw:
                data_dict[keyword] = [qstat_dict[target_job][keyword]]

        # Finish converting dynamic data into a pandas dataframe
        df = pd.DataFrame.from_dict(dict(data_dict))
        return df

    def to_csv(self, file, qstat_dict, target_job, overwrite=False):
        """
        Convert a qstat dictionary to a csv file that contains data, which can be
        plotted.

        :param file:  The absolute path of the csv file.
        :type file:  str.
        :param qstat_dict:  Qstat data that has been parsed into a dictionary.
        :type qstat_dict:  dict.
        :param target_job:  The target job that's being analyzed.
        :type target_job:  str.
        :param overwrite:  A flag to determine weather the csv file will be overwritten.
        :type overwrite:  bool.
        """
        data_file = Path(file)
        target_df = self.to_dataframe(qstat_dict=qstat_dict, target_job=target_job)
        if data_file.is_file():
            if not overwrite:
                with open(data_file, 'a') as _f:
                    target_df.to_csv(_f, header=False, index=False, index_label=False)
            else:
                with open(data_file, 'w') as _f:
                    target_df.to_csv(str(data_file), index=False, index_label=False)
        else:
            with open(data_file, 'w') as _f:
                target_df.to_csv(str(data_file), index=False, index_label=False)

    def to_sqlite(self):
        # Have a table that consists of static data per job,
        # and another table that keeps up with dynamic data.
        # Each dynamic row will have a primary key (the JobID)
        # that links to a static row.
        pass


class Qstat(BaseQstat):

    def __init__(self, job, wait_time=120, **kwargs):

        super().__init__(job=job, **kwargs)
        self.wait_time = wait_time
        self.watch_count = 0

    def countdown(self, wait_time=None):
        """
        This method takes a wait time and preforms a countdown in the terminal
        during the wait period.

        :param wait_time:  The time in seconds to wait.
        :type wait_time:  int.
        """
        while wait_time > 0:
            sys.stdout.write('Countdown: ' + '\033[91m' + str(wait_time) + '\033[0m' + '     \r')
            wait_time -= 1
            sleep(1)
        sys.stdout.write('\r')

    def watch(self, count=None, python_datetime=datetime.now(), max_count=None):
        """
        This watch method is directly used by the end user to collect data on the
        target job over time.

        :param count:  The count of qstat data points collected.
        :type count:  int.
        :param python_datetime:  A date time that will be added to the qstat data.
        :type python_datetime:  datetime.
        :param max_count:  The maximum number of times that Qstat will collect data
        points on the target job.
        :type max_count:  int.
        """
        if count is None:
            self.watch_count = 0
        self._watch(count=count, python_datetime=python_datetime, max_count=max_count)

    def _watch(self, count=None, python_datetime=datetime.now(), first_time=None, max_count=None):
        """
        This method should normally not be used by the end user.  It also collects
        data on the target job over time.

        :param count:  The count of qstat data points collected.
        :type count:  int.
        :param python_datetime:  A date time that will be added to the qstat data.
        :type python_datetime:  datetime.
        :param first_time:  A flag that determines weather it's the first time
        this function has been run for the target job.
        :type first_time:  bool.
        :param max_count:  The maximum number of times that Qstat will collect data
        points on the target job.
        :type max_count:  int.
        """
        # Count the number of data-points that have been taken during the watch.
        if count is None:
            self.watch_count += 1
        elif count is not None:
            self.watch_count = count + 1

        if first_time is None:
            first_time = True
        else:
            first_time = first_time
        try:

            self.run_qstat(csv_flag=True, sqlite_flag=False)
            self.qstat_log.info("Added data-point %s from qstat for %s." % (self.watch_count, self.target_job))
            if not first_time:
                if max_count is not None:
                    raise TargetJobKeyError
                self.countdown(wait_time=self.wait_time)
            self._watch(python_datetime=python_datetime, first_time=False, max_count=max_count)
        except TargetJobKeyError:
            if first_time:
                raise TargetJobKeyError("The target job cannot be found:  %s" % self.target_job)
            elif max_count is not None:
                self.qstat_log.info('Watched %s for %s iterations.' % (self.target_job, max_count))
            else:
                self.qstat_log.info('Finished watching %s' % self.target_job)

    def plot_data(self, data_file):
        """
        This method plots qstat data over time.  Three different interactive plotly
        line graphs are created.  One with mem and vmem, one with cpu percentage,
        and one with cpu time.

        :param data_file:  A csv file that contains qstat data.
        :type data_file:  str.
        """
        df = pd.DataFrame.from_csv(str(data_file))
        df_home = Path(data_file).parent

        dt = df["datetime"]
        job_state = list(df['job_state'])
        # Memory data
        ru_mem = list(df['resources_used.mem'])
        ru_mem = [int(str_num.replace("kb", "")) for str_num in ru_mem]
        ru_vmem = list(df['resources_used.vmem'])
        ru_vmem = [int(str_num.replace("kb", "")) for str_num in ru_vmem]
        # CPU data
        ru_cpupercent = list(df["resources_used.cpupercent"])
        ru_cpupercent = [int(str_num) for str_num in ru_cpupercent]
        ru_cput = list(df["resources_used.cput"])
        ru_cput = [int(str_num) for str_num in ru_cput]

        # Memory traces
        vmem_trace = go.Scatter(
            x=dt,
            y=ru_vmem,
            text=job_state,
            textposition='top center',
            mode='lines+text',
            name="Virtual Memory",
            line=dict(
                color='rgb(205, 12, 24)',
                width=4
            )
        )

        mem_trace = go.Scatter(
            x=dt,
            y=ru_mem,
            text=job_state,
            textposition='top center',
            mode='lines+text',
            name="Memory",
            line=dict(
                color='rgb(22, 96, 167)',
                width=4
            )
        )
        # CPU traces
        cpupercent_trace = go.Scatter(
            x=dt,
            y=ru_cpupercent,
            text=job_state,
            textposition='top center',
            mode='line+text',
            name="CPU Percentage",
            line=dict(
                color='rgb(205, 12, 24)',
                width=4
            )
        )

        cput_trace = go.Scatter(
            x=dt,
            y=ru_cput,
            text=job_state,
            textposition='top center',
            mode='line+text',
            name="CPU Time",
            line=dict(
                color='rgb(205, 12, 24)',
                width=4
            )
        )

        # Memory plot
        mem_html_file = df_home / "mem-plot.html"
        plotly.offline.plot([mem_trace, vmem_trace], filename=str(mem_html_file), auto_open=False)

        # CPU plots
        cpupercent_html_file = df_home / "cpupercent-plot.html"
        cput_html_file = df_home / "cput-plot.html"
        plotly.offline.plot([cpupercent_trace], filename=str(cpupercent_html_file), auto_open=False)
        plotly.offline.plot([cput_trace], filename=str(cput_html_file), auto_open=False)


class MultiQstat(object):

    def __init__(self, jobs, config_home="~/.pbs", cmd="qstat -f"):
        """
        This object is used to watch multiple pbs jobs and collect data on these
        jobs asynchronously.

        :param jobs:  A list of pbs jobs to watch.
        :type jobs:  list.
        :param config_home:  A configuration directory for storing pbs jobs.
        :type config_home:  str.
        :param cmd:  The qstat command used to produce the job status report.
        :type cmd:  str.
        """

        self.cmd = cmd
        self.target_jobs = jobs
        self.config_home = Path(config_home)
        self.job_dict = {}
        self.job_list = []

    def watch(self, jobs, infile=None, outfile=None, cmd=None, wait_time=120):
        """
        The watch method populates the job dictionary, which is used to watch
        the target jobs.  The most up to date job info is assigned to the
        job list class variable.

        :param jobs:  A list of target jobs.
        :type jobs:  list.
        :param infile:  The input file and the output file are used in tandem to determine the
        data file that will be used.  If only one of these values are given (infile/outfile), then
        it will be used as the data file.  If neither of these values are given, then a default file
        ("job_data.csv") will be used.  If both are given, then the infile data is appended to the outfile,
        which is used as the data file.
        :type infile:  str.
        :param outfile:  See infile.
        :type outfile: str.
        :param cmd:  The qstat command used to produce the job status report.
        :type cmd:  str.
        :param wait_time:  The amount of time to wait in between each point of data being collected.
        :type wait_time:  int.
        """
        self.job_dict = self.get_qstat_dict(jobs, infile=infile, outfile=outfile, cmd=cmd, wait_time=wait_time)
        self.job_list = self.multi_watch(job_dict=self.job_dict)

    def get_qstat_dict(self, jobs, infile=None, outfile=None, cmd=None, wait_time=120):
        """
        This job populates a dictionary with qstat objects that correspond to
        individual pbs jobs.

        :param jobs:  A list of target jobs.
        :type jobs:  list.
        :param infile:  The input file and the output file are used in tandem to determine the
        data file that will be used.  If only one of these values are given (infile/outfile), then
        it will be used as the data file.  If neither of these values are given, then a default file
        ("job_data.csv") will be used.  If both are given, then the infile data is appended to the outfile,
        which is used as the data file.
        :type infile:  str.
        :param outfile:  See infile.
        :type outfile: str.
        :param cmd:  The qstat command used to produce the job status report.
        :type cmd:  str.
        :param wait_time:  The amount of time to wait in between each point of data being collected.
        :type wait_time:  int.
        :return:  The job dictionary is returned.
        :rtype:  dict.
        """

        job_dict = {}
        for job in jobs:
            # Get qstat parameters for each target job
            home = str(self.config_home / job)
            if infile is None:
                infile = str(job) + '.csv'
            if outfile is None:
                outfile = str(job) + '.csv'
            _qstat = Qstat(job=job, home=home, infile=infile, outfile=outfile, cmd=cmd, wait_time=wait_time)
            # Create a dictionary value for each job
            job_dict[job] = _qstat
        return job_dict

    def multi_watch(self, job_dict):
        """
        The multi_watch method takes a dictionary of qstat jobs and
        uses them to asynchronously watch the target jobs in the dictionary.

        :param job_dict:  A dictionary whose values are Qstat objects.
        :type job_dict:  dict.
        :return:  A list of updated qstat jobs that are returned from the asynchronous
        watch method.
        :rtype:  list.
        """

        # Set up list/dict objects for saving qstat data
        tasks = []
        ioloop = asyncio.get_event_loop()

        for _qstat in job_dict.keys():
            # Append task list for asnychronous programming
            tasks.append(asyncio.ensure_future(self._async_watch(qstat=_qstat, count=_qstat.watch_count)))
        # Run task list and then close
        job_list = ioloop.run_until_complete(asyncio.wait(tasks))
        ioloop.close()
        return job_list

    async def _async_watch(self, qstat, first_time=None, count=None):
        """
        This asynchronous method runs Qstat commands for multiple jobs.  The
        asynchronous component is used during the waiting period in between data
        points.

        :param qstat:  A called qstat object.
        :type qstat:  Qstat.
        :param first_time:  A flag that helps determine if it's the first time that a function is used.
        :type first_time:  bool.
        :param count:  A number that corresponds to the number of times that a job has been observed.
        :type count:  int.
        :return:  An updated Qstat object.
        :rtype:  Qstat.
        """
        # Count the number of data-points that have been taken during the watch.
        if count is None:
            qstat.watch_count += 1
        elif count is not None:
            qstat.watch_count = count + 1

        if first_time is None:
            first_time = True
        else:
            first_time = first_time

        try:
            qstat.run_qstat(csv_flag=True, sqlite_flag=False)
            qstat.qstat_log.info("Added data-point %s from qstat for %s." % (qstat.watch_count, qstat.target_job))
            if not first_time:
                await asyncio.sleep(qstat.wait_time)
            temp_qstat = self._async_watch(qstat=qstat, first_time=False)
        except TargetJobKeyError:
            if first_time:
                raise TargetJobKeyError("The target job cannot be found:  %s" % qstat.target_job)
            else:
                qstat.qstat_log.info('Finished watching %s' % qstat.target_job)
                temp_qstat = qstat
        qstat = temp_qstat
        return qstat

