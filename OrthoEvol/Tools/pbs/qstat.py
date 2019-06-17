import os
import shutil
import click
import subprocess
import yaml
from pathlib import Path
import pandas as pd
from collections import OrderedDict
from dateutil import parser
import random
from datetime import datetime
from time import sleep
from OrthoEvol.Manager.config import scripts, yml
from pkg_resources import resource_filename
from OrthoEvol.Tools.pbs import _setup_yaml
from OrthoEvol.Tools import pbs


class BaseQwatch(object):
    # Setup the yaml with OrderedDict
    _setup_yaml()
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

    def __init__(self, jobs=None, users=os.getlogin(), email=None, infile=None, watch=None, filename_pattern=None,
                 directory=f"qwatch{random.randint(9001, 100000)}", cmd="qstat -f", sleeper=120,
                 clean=False, slack=None):

        # TODO-ROB: Implement emial notifications
        # TODO-ROB: Implement slack notifications
        # TODO-ROB: Add logger (through click?)
        # TODO-ROB: Do the rework in get_dicts()
        """
        A class for parsing "qstat -f" output on SGE systems for monitoring
        jobs and for making smart downstream decisions about resource allocation.
        This class uses asyncio in order to monitor multiple jobs at the same time.

        :param jobs:  Jobs to look for.
        :type jobs: list
        :param email: Email to send files to
        :type email: list
        :param infile: A file that contains 'qstat -f' output
        :type infile: str
        :param watch: A flag that indicates that the selected jobs are to be watched
        :type watch: bool
        :param filename_pattern: A pattern used to create all of the files.
        :type filename_pattern: str
        :param directory: The directory that is used to contain the output files.
        :type directory: str
        :param users: Users to look for in the job que.
        :type users: list
        :param cmd: The qstat command to use.  This should only be able to change if new parsers are created.
        :type cmd: str
        :param sleeper: The minimum length of time in between each data point.
        :type sleeper: int
        """
        self.cmd = cmd
        # Get a user list
        if not users:
            self.users = []
        elif isinstance(users, str):
            self.users = [users]
        elif isinstance(users, list):
            self.users = users
        elif isinstance(users, tuple):
            self.users = list(users)
        else:
            raise TypeError("The users parameter is a single user string or a multi-user list.")
        # Get a job list
        if not jobs:
            self.jobs = []
        elif isinstance(jobs, str):
            self.jobs = [jobs]
        elif isinstance(jobs, list):
            self.jobs = jobs
        elif isinstance(jobs, tuple):
            self.users = list(jobs)
        else:
            raise TypeError("The jobs parameter is a single job string or a multi-job list.")

        # Initialize file types
        self.qstat_filename = Path()
        self.yaml_filename = Path()
        self.data_filename = Path()
        self.info_filename = Path()
        self.plot_filename = Path()
        self._new_keyword_log = Path()
        # Initialize utility files
        self._yaml_config = resource_filename(yml.__name__, 'qstat_dict.yml')
        self.r_line_graph_filename = resource_filename(scripts.__name__, 'line_graph_workflow.R')
        self._async_filename = resource_filename(pbs.__name__, 'watcher.py')
        # Initialize other file based attributes
        self.infile = infile
        while directory.is_dir():
            directory = Path(f"qwatch{random.randint(9001, 100000)}")
        self.directory = directory
        self.directory.mkdir(parents=True)
        self._temp_yml = directory / Path("temp_yml.yml")
        self.filename_pattern = filename_pattern
        self.initialize_file_names()

        # Initialize other attributes
        self.orig_jobs = self.jobs
        self.email = email
        self.watch = watch
        self.sleeper = sleeper

    def initialize_file_names(self):
        # Get a filename pattern based on other user input
        if not self.filename_pattern:
            if self.infile:
                filename_pattern = f"{self.infile}"
            elif len(self.users) == 1 and len(self.jobs) == 0:
                filename_pattern = f"{self.users[0]}"
            elif len(self.jobs) == 1 and len(self.users) == 0:
                filename_pattern = f"{self.jobs[0]}"
            else:
                current_user = os.getlogin()
                _id = random.randint(10000, 99999)
                filename_pattern = f"{current_user}_{_id}"
            self.filename_pattern = filename_pattern.replace('.', '_')

        # Create file names using the pattern
        self.yaml_filename = Path(self.directory) / Path(f"{self.filename_pattern}.yml")
        self.data_filename = Path(self.directory) / Path(f"{self.filename_pattern}.data")
        self.info_filename = Path(self.directory) / Path(f"{self.filename_pattern}.info")
        self.plot_filename = Path(self.directory) / Path(f"{self.filename_pattern}_plot.png")
        self._new_keyword_log = self.directory / Path('lost_and_found_qstat_keywords.log')

    def save_qstat_data(self):
        """
        Save the 'qstat -f' output to the qstat_file or set the infile to the qstat_file.
        """
        if self.infile:
            self.qstat_filename = Path(self.infile)
        else:
            self.qstat_filename = Path(self.directory) / Path(f"{self.filename_pattern}.qstat")
            qstat = subprocess.Popen(self.cmd, stderr=subprocess.PIPE,
                                     stdout=subprocess.PIPE, shell=True,
                                     encoding='utf-8', universal_newlines=False)
            out = qstat.stdout.readlines()
            error = qstat.stderr.readlines()
            with open(self.qstat_filename, 'w') as qf:
                qf.writelines(out)
            print(error)

    def qstat_to_filtered_yaml(self):
        """
        This function saves the qstat information to a yaml file.  It parses the qstat data, filters out the unwanted
        jobs, and then validates the data before saving it in YAML format.  It doesn't return anything, but it does
        create a file.
        """

        job_dict = self.qstat_to_dict()

        # Filter and keep only the selected jobs and then create a YAML file
        kept_jobs, kept_dict = self.filter_jobs(job_dict)
        self.jobs = kept_jobs

        # If the yaml file doesnt exist then update the jobs and dump the qstat data
        if not self.yaml_filename.is_file() or self.watch is True:
            with open(self.yaml_filename, 'w') as yf:
                yaml.dump(kept_dict, stream=yf, default_flow_style=False)

        # If the yaml file is empty, then overwrite it.
        else:
            with open(self.yaml_filename, 'r') as yf:
                test = yaml.load(yf)
                if not test:
                    with open(self.yaml_filename, 'w') as yf2:
                        yaml.dump(kept_dict, stream=yf2, default_flow_style=False)

    def filtered_yaml_to_dict(self):
        """
        This function loads the yaml file created with qstat_to_filtered_yaml and validates the qstat data.
        :return: A dictionary of jobs.
        :rtype: dict
        """
        if not self.yaml_filename.is_file() or self.watch is True:
            self.qstat_to_filtered_yaml()
        else:
            with open(self.yaml_filename, 'r') as yf:
                test = yaml.load(yf)
                if not test:
                    self.qstat_to_filtered_yaml()
        with open(self.yaml_filename, 'r') as yf2:
            jobs_dict = yaml.load(yf2)
        return jobs_dict

    def update_qstat_data(self, save, process, data):
        # Save qstat output
        if save:
            self.save_qstat_data()
        # Process qstat output and create the filtered yaml file
        if process:
            self.qstat_to_filtered_yaml()
        # Take the filtered yaml file and create a dictionary
        if data:
            _data = self.filtered_yaml_to_dict()
            return _data

    def qstat_to_dict(self):
        """
        The qstat parser takes the qstat file from the user infile or from the
        'qstat -f' command and parses it.  It uses the qstat keywords found in the
        qstat yaml file.
        :return:  A dictionary of jobs.
        :rtype: dict
        """
        mast_dict = OrderedDict()
        job_count = 0
        phrase_continuation_flag = None
        with open(self._yaml_config, 'r') as yf:
            qstat_keywords = yaml.load(yf)
        with open(self.qstat_filename, 'r') as qf:
            qstat_sentence = None
            continuation_phrase = ""
            qstat_phrase = ""
            prev_item = None
            for item in qf.readlines():
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

    def get_dicts(self, python_datetime=None):
        """
        This function takes the filtered yaml file and creates a nested dictionary object that is suitable for
        creating/updating the info and data csv files.

        :param python_datetime:
        :type python_datetime:
        :return:
        :rtype:
        """
        df = self.filtered_yaml_to_dict()
        # jobs_dict = self.filtered_yaml_to_dict()
        #df = OrderedDict()
        master_dict = OrderedDict()
        info_dict = OrderedDict()
        data_dict = OrderedDict()
        # TODO-ROB: Rework this or rework filterd_yaml_to_dict
        # for job in jobs_dict.keys():
        #     row = OrderedDict()
        #     _job = OrderedDict()
        #     row = jobs_dict[job]
        #     _job[job] = row
        #     df.update(_job)
        for job in df.keys():
            var_dict = OrderedDict()
            # Store PBS Environment Variables
            for var in df[job]['Variable_List'].keys():
                var_dict[var] = [df[job]['Variable_List'][var]]
            info_dict[job] = var_dict
            data_dict[job] = OrderedDict()
            # Store the python datetime
            if python_datetime:
                info_dict[job]["datetime"] = [python_datetime]
                data_dict[job]["datetime"] = [python_datetime]
            # Loop through the keywords
            for keyword in df[job].keys():
                # Store all of the dynamic data
                if keyword in self.__static_kw:
                    if keyword != "Variable_List":
                        # Store the job time values in a datetime format
                        if keyword in self.__job_time_kw:
                            info_dict[job][keyword] = [str(parser.parse(df[job][keyword]))]
                        else:
                            info_dict[job][keyword] = [df[job][keyword]]
                # Store all of the dynamic data
                elif keyword in self.__dynamic_kw:
                    data_dict[job][keyword] = [df[job][keyword]]

        master_dict["info"] = info_dict
        master_dict["data"] = data_dict
        return master_dict

    def get_dataframes(self, python_datetime=None):
        """
        This function returns a dataframe of the current yaml data.  It is specifically used by the update_csv function
        for creating csv files of the qstat data
        :param python_datetime:
        :type python_datetime:
        :return:
        :rtype:
        """
        master_dict = self.get_dicts(python_datetime=python_datetime)
        master_df = OrderedDict()
        for key in master_dict.keys():
            if len(self.jobs) > 1:
                master_df[key] = pd.DataFrame.from_dict(master_dict[key])
            else:
                for job in master_dict[key].keys():
                    master_df[key] = pd.DataFrame.from_dict(dict(master_dict[key][job]))

        return master_df

    def get_info(self, data_frame=False, python_datetime=None):
        """
        Get the data in the info file as a python object (dictionary or dataframe).
        :param data_frame:
        :type data_frame:
        :param python_datetime:
        :type python_datetime:
        :return:
        :rtype:
        """
        if data_frame:
            _data = self.get_dataframes(python_datetime=python_datetime)
        else:
            _data = self.get_dicts(python_datetime=python_datetime)
        return _data["info"]

    def get_data(self, data_frame=False, python_datetime=None):
        """
        Get the data in the data file as a python object (dictionary or dataframe).
        :param data_frame:
        :type data_frame:
        :param python_datetime:
        :type python_datetime:
        :return:
        :rtype:
        """
        if data_frame:
            _data = self.get_dataframes(python_datetime=python_datetime)
        else:
            _data = self.get_dicts(python_datetime=python_datetime)
        return _data["data"]

    def filter_jobs(self, job_dict):
        """
        Filter a job dict based on the input user or job list.
        :param job_dict:
        :type job_dict:
        :return:
        :rtype:
        """
        kept_jobs = []
        if len(self.jobs) != 0 or len(self.users) != 0:
            # Create a new list of kept jobs
            for j in job_dict.keys():
                if len(self.users) > 0:
                    if job_dict[j]["Variable_List"]["PBS_O_LOGNAME"] in self.users:
                        kept_jobs.append(j)
                elif len(self.jobs) > 0:
                    if job_dict[j]["Job_Id"] in self.jobs:
                        kept_jobs.append(j)
        # If no user or job is given then use all the jobs
        else:
            kept_jobs += list(job_dict.keys())
        # Format the kept job data
        kept_jobs = list(set(kept_jobs))
        kept_dict = OrderedDict((k, job_dict[k]) for k in kept_jobs)
        return kept_jobs, kept_dict

    def clean_output(self, ext_list=[".info", ".data", ".png", ".log"]):

        for root, dirs, files in os.walk(self.directory):
            for file in files:
                if Path(file).suffix not in ext_list:
                    file_to_delete = str(Path(root) / Path(file))
                    os.remove(file_to_delete)


class Qwatch(BaseQwatch):

    def __init__(self, **kwargs):
        """
        This is the top level qwatch class.  It watches jobs, updates csv data files, and plots
        graphs by calling an R script file.  If multiple jobs are being watched, then asyncio
        is utilized for efficiency.
        :param kwargs:
        :type kwargs:
        """
        super().__init__(**kwargs)
        self.qwatch = BaseQwatch

    def _get_subset_kwargs(self, skipped_kwargs):
        """
        Get kwargs that are subset from the original object.  This is used specifically for the async_qwatch.py
        script.
        :param skipped_kwargs:
        :type skipped_kwargs:
        :return:
        :rtype:
        """
        _kwargs = {}
        for var, attr in self.__dict__.items():
            if var not in skipped_kwargs:
                _kwargs[var] = attr
        return _kwargs

    def new_keyword_logger(self, data):
        """
        This logger function checks for new keywords, and adds them to a
        lost_and_found_keywords.log file.
        :param data: The data to check.
        """
        with open(self._yaml_config, 'r') as yc:
            _checker_dict = yaml.load(yc)
            _checker_list = list()
            _checker_list.append('Job Id')
            _checker_list = _checker_list + list(_checker_dict['Job Id'].keys())
            _checker_list = _checker_list + list(_checker_dict['Job Id']['Variable_List'].keys())
        _data_index_list = list(data.index)
        diff = list(set(_data_index_list) - set(_checker_list))
        if diff == [0]:
            diff = []
        if len(diff) != 0:
            if not self._new_keyword_log.is_file():
                with open(self._new_keyword_log, 'w') as cl:
                    cl.write("***********************************   NOTICE  ***********************************\n\n")
                    cl.write("This is a log file for new keywords that were missed by the parser.")
                    cl.write("Please email this file to Rob Gilmore at rgilmore@umc.edu, if new keywords were logged.")
                    cl.write("Thanks!\n\n")
                    cl.write("***********************************   NOTICE  ***********************************")
            with open(self._new_keyword_log, 'a') as cl:
                cl.write('Unresolved qstat keyword: %s' % diff)

    def update_csv(self, file, data):
        """
        The update_csv file is named appropriately, but it also checks for undocumented qstat keywords.
        :param file:
        :type file:
        :param data:
        :type data:
        :return:
        :rtype:
        """

        # Check for undocumented keywords
        self.new_keyword_logger(data)

        # Append new data to the csv file (either .info or .data
        if file.is_file():
            file_df = pd.read_csv(file, index_col=False)
            updated_df = pd.concat([file_df, data])
            updated_df.to_csv(file, index=False, index_label=False)
        else:
            data.to_csv(file, index=False, index_label=False)

    def watch_jobs(self):
        """
        This function watches all of the jobs that passed the filter.  Multiple jobs will
        utilize asyncio and the seperate async_qwatch.py script.
        :return:
        :rtype:
        """
        # Initialize data
        self.save_qstat_data()
        self.qstat_to_filtered_yaml()
        if len(self.jobs) > 1:
            kw_dict = {}
            kw_dict["jobs"] = self.jobs
            kw_dict["kwargs"] = self._get_subset_kwargs(skipped_kwargs=["jobs", "directory", "qwatch", "watch", "_yaml_config",
                                                                        "info_filename", "plot_filename", "qstat_filename",
                                                                        "data_filename", "yaml_filename", "users",
                                                                        "filename_pattern", "orig_jobs",
                                                                        "_async_filename", "r_line_graph_filename"])

            kw_dict["sleeper"] = 5
            kw_dict["directory"] = self.directory
            # Call the asyncio setup for multiple jobs and begin watching
            with open(self._temp_yml, 'w') as ty:
                yaml.dump(kw_dict, stream=ty, default_flow_style=False)

            qstat = subprocess.Popen(f'python3.6 {self._async_filename} -yamlfile {self._temp_yml}', stderr=subprocess.PIPE,
                                     stdout=subprocess.PIPE, shell=True,
                                     encoding='utf-8', universal_newlines=False)
            out = qstat.stdout.read()
            error = qstat.stderr.read()
            print(f'out: {out}')
            print(f'error: {error}')
        elif len(self.jobs) == 1:
            # Call the single job watcher
            self._watch(datetime.now())

    def _watch(self, python_datetime=None, first_time=True):
        """Wait until a job finishes and get updates."""

        job_dict = self.update_qstat_data(save=True, process=True, data=True)

        if job_dict:
            data = self.get_data(data_frame=True, python_datetime=python_datetime)
            self.update_csv(file=self.data_filename, data=data)
            print(f"Updated qstat data for {job_id}")

            if job_dict[self.jobs[0]]['job_state'] == 'Q':
                sleep(self.sleeper)
                self._watch(datetime.now(), first_time=True)
            elif job_dict[self.jobs[0]]['job_state'] == 'R':
                # Create the static data file on the first instance of a running job
                if first_time:
                    info = self.get_info(data_frame=True, python_datetime=python_datetime)
                    self.update_csv(file=self.info_filename, data=info)
                sleep(self.sleeper)
                self._watch(datetime.now(), first_time=False)

        return f'Finished {self.jobs[0]}'

    def plot_memory(self, data_file=None, info_file=None, file_pattern=None, directory=None, jobs=None, rdata_save=False):
        """
        Plot a graph of the memory vs time.  Utilizes an R script which creates a ggplot2 line graph.
        :param data_file:
        :type data_file:
        :param info_file:
        :type info_file:
        :param file_pattern:
        :type file_pattern:
        :param directory:
        :type directory:
        :param jobs:
        :type jobs:
        :param rdata_save:
        :type rdata_save:
        :return:
        :rtype:
        """
        if jobs:
            for job in jobs:
                print(f'job: {job}')
                if directory:
                    dir_path = Path(directory) / Path(job)
                else:
                    dir_path = Path(job)
                print(f'directory: {dir_path}')
                if not data_file:
                    if not file_pattern:
                        df = dir_path / Path(f'{job}.csv')
                    else:
                        df = dir_path / Path(f'{file_pattern}.csv')
                print(f'data_file: {df}')
                if not info_file:
                    if not file_pattern:
                        inf = dir_path / Path(f'{job}_info.txt')
                    else:
                        inf = dir_path / Path(f'{file_pattern}.txt')
                print(f'info_file: {inf}')
                cmd = f'Rscript {self.r_line_graph_filename} -d {str(df)} -i {str(inf)}'
                if rdata_save:
                    cmd = cmd + ' --rdata_save'
                if file_pattern:
                    cmd = cmd + f' --name {file_pattern}'
                print(f'cmd: {cmd}')
                plot = subprocess.Popen(cmd,
                                        stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True,encoding='utf-8',
                                        universal_newlines=False)
                out = plot.stdout.readlines()
                error = plot.stderr.readlines()
                print(out)
                print(error)
        else:
            plot = subprocess.Popen(f'Rscript {self.r_line_graph_filename} -d {str(self.data_filename)} '
                                    f'-i {str(self.info_filename)} --name {self.filename_pattern}', stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE, shell=True, encoding='utf-8', universal_newlines=False)
            out = plot.stdout.readlines()
            error = plot.stderr.readlines()
            print(out)
            print(error)
