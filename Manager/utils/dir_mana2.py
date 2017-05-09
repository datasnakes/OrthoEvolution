##############################################################################
# PyCharm Community Edition
# -*- coding: utf-8 -*-
"""
Orthologs-Project
Directory_management updated on 11/15/2016 at 11:30 AM
##############################################################################

    Input:  A string that represents a path for the main file systems

    Output:  Custom class variables for import in the project file template.

    Description:  This file manges the directories for projects created in PyCharm.
    Each new project gets it's own function.

##############################################################################
@author: rgilmore
"""
##############################################################################
# Libraries:


import json
import os
import shutil
import time
from pathlib import Path
from cookiecutter.main import cookiecutter
from Manager.utils import treelib2
from Manager.utils.json_to_newick import _parse_json
import Cookies
import Manager
import Orthologs
import Tools

#from project_mana import project_mana  # //TODO-ROB: Add a configure function to proj_mana to get the root directory using the project name
##############################################################################
# Directory Initializations:

class dir_mana(object):
    """This class organizes a directory tree for a project.
    The dir_mana() class will help with organization and it will
    help to instantly access the proper directories on command.

    It is advised to set up a function that is named after the project.
    Each new project will be called from the project() function.

    See GPCR-Orthologs-Project."""

# //TODO-ROB Add JSON loading of the different directory variables
# //TODO-Rob change project to projects and add another variable called project_type

    def __init__(self, home=os.getcwd(), proj_mana="proj_mana.json"):
        """Initialize the directory tree for the project.
        Each project will have a home directory in addition to the following:
        ."""
        # config = tablib.Dataset().load(open('config.yaml').read())
        self.__file_home = Path(home)  # Home of the file calling this class
        # TODO-ROB The user information will only be accessible via the flask user model.
        # #TODO-ROB This is a helper class for FLASK
        # Below are the PyPi path strings
        #    The first group is to access the cookiecutter templates
        self._Cookies = Path(Cookies.__path__.path[0])
        self._new_repo = self._Cookies / Path('new_repository')
        self._new_user = self._Cookies / Path('new_user')
        self._new_project = self._Cookies / Path('new_project')
        self._new_research = self._Cookies / Path('new_research')
        self._new_app = self._Cookies / Path('new_app')
        #    The second group is for the Manager module
        self._Manager = Path(Manager.__path__.path[0])
        self._index = self._Manager / Path('index')
        self._utils = self._Manager / Path('utils')
        self._shiny = self._Manager / Path('shiny')
        #    The third group is for the Orthologs module
        self._Orthologs = Path(Orthologs.__path__.path[0])
        # TODO-ROB Add the other paths here
        self._Tools = Path(Tools.__path__.path[0])

    def flask_user_config(self):
        if

    # Map the main project directory.
    def get_dir_map(self, top, ignore=None):
        default_ignore = ['.git', '.idea']
        if ignore is not None:
            ignore += default_ignore
        else:
            ignore = default_ignore
        # Treelib will help to map everything and create a json at the same time
        tree = treelib2.Tree()
        tree.create_node('.', top)
        # Walk through the directory of choice (top)
        # Topdown is true so that directories can be modified in place
        for root, dirs, files in os.walk(top, topdown=True):
            # Only remove directories from the top
            if root == top:
                print(root)
                try:
                    dirs[:] = set(dirs) - set(ignore)  # Remove directories from os.walk()
                    print(dirs)
                except ValueError:
                    pass
            for d in dirs:
                rd = str(Path(root) / Path(d))
                tree.create_node(d, identifier=rd, parent=root)
            for f in files:
                tree.create_node(f, parent=root)
        return tree

    # def user_dir_config(self, username):
    #     # TODO-ROB USE SQL here to see if the user db contains the username
    #     if username in os.listdir(self.users):
    #         return UserWarning('This user already exists')
    #     else:
    #         # TODO-ROB CREATE THESE IN A VIRTUAL ENVIRONMENT FOR EACH USER
    #         # TODO-ROB The virtual environment can be the name of the user
    #         # TODO-ROB When the user logs in, they will activate the virtual environment
    #         for path in self.user_dict.values():
    #             Path.mkdir(path, parents=True, exist_ok=True)
    #         for path in self.user_project_dict.values():
    #             Path.mkdir(path, parents=True, exist_ok=True)
    #
    # def user_project_config(self, project_name, private=False):
    #     if private is False:
    #         path = self.user_project_dict['public']
    #     else:
    #         path = self.user_project_dict['private']
    #     Path.mkdir(path / Path(project_name))

    def get_newick_dir_map(self, top, ignore=None):
        """Takes a treelib tree created by get_dir_map and returns
        a tree a dir_map in newick format.  This will be useful for Bio.Phylo
        applications."""

        tree = self.get_dir_map(top, ignore)
        Ntree = _parse_json(tree.to_jsonnewick())
        return Ntree

    def proj_mana(self, js_file):
        with open(js_file, 'r') as file:
            pdr = json.load(file)
        return pdr

    # # //TODO-ROB Find a different way to return a
    # def path_list_make(self, path, o_path=None):
    #     # Takes a path and reduces it to a list of directories within the project
    #     # An optional attribute (o_path) is give so that a deeper path within the project can be used
    #     home = str(self.__project_home).split('/')
    #     path_list = str(path).split('/')
    #     for item in home:
    #         if item in path_list:
    #             path_list.remove(item)
    #     # path_list = set(p) - set(home)
    #     if o_path is not None:
    #         o_path = str(o_path).split('/')
    #         for item in o_path:
    #             if item in path_list:
    #                 path_list.remove(item)
    #         # path_list = set(path_list) - set(o_path)
    #     return path_list

    # //TODO-ROB utilize Path.mkdir(parents=TRUE) instead
    # def dir_make(self, path, path_list):
    #     # Takes a path list which is a list of folder names
    #     # path_list created by str(path).split('/')
    #     # The path_list appends to path, which is already an established directory inside the project
    #     t = None
    #     for item in path_list:
    #
    #         if os.path.isdir(path + '/' + item): # If for some reason the directory already exists...
    #             path += '/' + item  # Append a directory
    #             continue
    #         path += '/' + item  # Append a directory
    #         os.mkdir(path)
    #     if len(os.listdir(path)) > 0:
    #         path, t = self.dir_archive(path, path_list='')
    #     return path, t
    #
    # # //TODO-ROB Change to using a compression module https://pymotw.com/2/compression.html
    # def dir_archive(self, path, path_list):
    #     # Use the path that you want to update/add to
    #     # Returns path and the time stamp (could be None)
    #     unique_dir = False
    #     archive_path = path
    #     for item in path_list:
    #
    #         path += '/' + item  # Append a directory
    #         if os.path.isdir(path):  # If the child directory exists
    #             archive_path = path  # Then update the dir_archive path and continue
    #             continue
    #         else:                     # If the child directory doesnt exist
    #             unique_dir = True     # Then raise the flag
    #             os.mkdir(path)        # And make a directory
    #
    #     if unique_dir is False:  # Only dir_archive if the final child directory is not unique (via unique_dir = False)
    #         t = time.strftime("%m%d%Y-%I%M%S")
    #         new_archive = self.Archive + '/' + t  # Creates a time stamped directory
    #
    #         os.mkdir(new_archive)
    #         for item in os.listdir(archive_path):
    #             if os.path.isfile(archive_path + '/' + item):  # Only dir_archive the FILES
    #                 shutil.move(archive_path + '/' + item, new_archive)
    #         return path, t
    #     else:
    #         return path, None





# TODO-ROB Add a instance that stores new paths inside of a text file and has an updating/overwriting ability

# TODO-ROB CALL ON THIS INSTNACE IN __INIT__ SO THAT VARIABLES CAN BE CREATED BASED ON THE PATHS









