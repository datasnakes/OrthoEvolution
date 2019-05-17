#!%s/python

from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Manager.database_dispatcher import DatabaseDispatcher
from OrthoEvol.utilities import FullUtilities
import yaml
import os
from time import sleep
from pathlib import Path
def sleeper(secs):
   sleep(secs)

def _dispatch_config(config_file):
   # Initialize variables
   utils = FullUtilities()
   nrr_dispatcher = {"NCBI_refseq_release": {"upload": []}}
   nrr_config = {"NCBI_refseq_release": {"upload": []}}
   
   # Templated variables
   file_list = %s
   pbs_dict = %s
   database_path = "%s"
   upload_number = %s
   email = "%s"
   
   # Config file variables
   with open(config_file, 'r') as cf:
      cd = yaml.load(cf)
      cd = cd['Database_config']['Full']['NCBI']['NCBI_refseq_release']
   seqformat = cd['seqformat']
   seqtype = cd['seqtype']
   collection_subset = cd['collection_subset']
   activate = cd['activate']

   # Get a list of files to upload
   file_list = os.listdir(database_path)
   file_list = [x for x in file_list if x.endswith(str(seqformat))]
   file_dict = {}
   for file in file_list:
      file_dict[file] = Path(database_path, file).stat().st_size
   sub_upload_lists = utils.group_files_by_size(file_dict, groups=upload_number)
   add_to_default = 0

   # Create configuration for dispatching PBS jobs
   for sub_list in sub_upload_lists:
      add_to_default += 1
      nrr_dispatcher["NCBI_refseq_release"]["upload"].append(utils.refseq_jobber)
      nrr_dispatcher["NCBI_refseq_release"]["upload"].append(sleeper)
      code_dict_string = str({
         "collection_subset": collection_subset,
         "seqtype": seqtype,
         "seqformat": seqformat,
         "upload_list": sub_list,
         "add_to_default": add_to_default
      })
      # Create a Python script for this in the package
      sge_code_string = \
         "from OrthoEvol.Manager.management import ProjectManagement\n" \
         "from OrthoEvol.Manager.database_dispatcher import DatabaseDispatcher\n" \
         "from OrthoEvol.Manager.config import yml\n" \
         "from pkg_resources import resource_filename\n" \
         "import yaml\n" \
         "pm_config_file = \"%%s\"\n" \
         "with open(pm_config_file, \'r\') as f:\n" \
         "   pm_config = yaml.load(f)\n" \
         "pm = ProjectManagement(**pm_config[\"Management_config\"])\n" \
         "code_dict_string = %%s\n" \
         "R_R = DatabaseDispatcher(config_file=\"%%s\", proj_mana=pm, **code_dict_string)\n" \
         "R_R.upload_refseq_release_files(collection_subset=\"%%s\", seqtype=\"%%s\", seqformat=\"%%s\", upload_list=%%s, add_to_default=%%s)\n" %% \
         (config_file, code_dict_string, config_file, collection_subset, seqtype, seqformat, sub_list, add_to_default)
      nrr_config["NCBI_refseq_release"]["upload"].append({
         "code": sge_code_string,
         "base_jobname": "upload_rr_%%s",
         "email_address": email,
         "id": add_to_default,
         "activate": activate,
         "config_dict": pbs_dict})
      nrr_config["NCBI_refseq_release"]["upload"].append({'secs':15})
   return nrr_dispatcher, nrr_config

# Setup project management and function dispatcher
config_file = "%s"
with open(config_file, 'r') as f:
   pm_config = yaml.load(f)
pm = ProjectManagement(**pm_config["Management_config"])
dd = DatabaseDispatcher(config_file, pm)

# Dispatch PBS jobs
disp, conf = _dispatch_config(config_file)
dd.dispatch(strategies=list(disp.keys()),dispatcher=disp, configuration=conf)