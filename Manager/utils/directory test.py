'''This file is to simulate what happens when a user logs in.
Environment variables are read.'''

import os
from pathlib import Path
home = os.getcwd()
users = 'Rob'
__file_home = Path(home)  # Home of the file calling this class

current_user = os.environ['ACTIVE_USER']  # TODO-ROB:  FLASK detail
user_dict = {}
user_path = users / Path(current_user)
###  THIS CONFIG CAN BE DONE UNDER ~/users/config.ini
"""Configuration for the users directory.  Each user gets this."""
user_dict['username'] = user_path
user_dict['index'] = user_path / Path('index')
user_dict['log'] = user_path / Path('log')
user_dict['manuscripts'] = user_path / Path('manuscripts')
user_dict['other'] = user_path / Path('other')
user_dict['projects'] = user_path / Path('projects')
"""Configuration for each users projects directory."""
project_status = 'other'
user_dict['public'] = user_dict['projects'] / Path('public')
user_dict['private'] = user_dict['projects'] / Path('private')
# TODO-ROB: Make a project category that allows users to make some research targets private and othes public
user_dict['other'] = user_dict['projects'] / Path('other')

## THIS CONFIG CAN BE DONE UNDER ~/users/$CURRENT_USER/projects/config.ini
current_dataset = 'GPCR'  # TODO-ROB os.environ['CURRENT_DATASET']
user_dict['data_sets'] = {}
user_dict['data_sets'][current_dataset] = user_dict[project_status] / Path(current_dataset)

current_research_target = 'comparative_genetics'  # TODO-ROB: os.environ['CURRENT_RT']

project_research_target_dict['']
if project_status is 'other':
    # RT_status = {}
    # RT_status['comparative_genetics'] = 'private'
    # RT_status['polymorphisms'] = 'public'
    RT_status = 'private'
    user_dict['data_sets'][current_dataset][RT_status] =
    user_dict['data_sets'][current_dataset] / Path(current_research_target)


dataset_dict = {}
dataset_dict['data_set'] = user_project_dict[project_status] / Path(current_dataset)
dataset_dict[]
