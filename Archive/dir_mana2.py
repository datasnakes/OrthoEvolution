import json
import os
import shutil
import time
from pathlib import Path
from datasnakes.Orthologs import treelib2
from Manager.dir_map.json_to_newick import _parse_json

class dir_mana(object):
    """This class organizes a directory tree for a project.
    The Mana() class will help with organization and it will
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

        self.current_user = os.environ['ACTIVE_USER']  # TODO-ROB:  FLASK detail
        self.user_dict = {}
        self.user_path = self.users / Path(self.current_user)
        self.user_dict['username'] = self.user_path
        self.user_dict['index'] = self.user_path / Path('index')
        self.user_dict['log'] = self.user_path / Path('log')
        self.user_dict['manuscripts'] = self.user_path / Path('manuscripts')
        self.user_dict['other'] = self.user_path / Path('other')
        self.user_dict['projects'] = self.user_path / Path('projects')

        self.project_status = ''
        self.user_project_dict = {}
        self.user_project_dict['public'] = self.user_dict['projects'] / Path('public')
        self.user_project_dict['private'] = self.user_dict['projects'] / Path('private')
        # TODO-ROB: Make a project category that allows users to make some research targets private and othes public
        self.user_project_dict['other'] = self.user_dict['projects'] / Path('other')
        self.current_dataset = ''
        self.user_project_dict[self.current_dataset] = self.user_project_dict[self.project_status] / Path(self.current_dataset)
        self.project_research_target_dict = {}
        if self.project_status is 'other':
            self.user_project_dict[self.current_dataset]


        self.dataset_dict = {}
        self.dataset_dict['data_set'] = self.user_project_dict[self.project_status] / Path(self.current_dataset)
        self.dataset_dict[]

        self.re
        #self.__venv_home = os.environ['VIRTUAL_ENV']
        #self.__python = self.__venv_home / Path('bin')

        #self.__project_management = self.proj_mana(proj_mana)  # dict {[projects] : {[data-sets] : [research targets]}}
        # if pdr is False:
        #     self.__project = project  # Project Name
        #     self.__project_home = Path()  # Initialize home to a PathLike object.
        #     self.__dataset = self.__project_home / Path(dataset)
        #
        #     # The project name is also the name of the project's main directory.
        #     # Make sure that the file's home is in the project's directory tree.
        #     if self.__project not in list(self.__file_home.parts):
        #         raise NotADirectoryError
        #     # If the file isn't in the directory tree throw an Error.
        #
        #     # Append the project_home variable with child directories until the
        #     # project's main directory is reached.
        #     for item in list(self.__file_home.parts):
        #         self.__project_home = self.__project_home / Path(item)
        #         if item == self.__project:
        #             break
        #     self.project()
        # else:
        #
        #     self.top_level_struct()
        #     self.project()
        #
        # if exists is False:
        #     self.create_dirs()

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

    def user_dir_config(self, username):
        # TODO-ROB USE SQL here to see if the user db contains the username
        if username in os.listdir(self.users):
            return UserWarning('This user already exists')
        else:
            # TODO-ROB CREATE THESE IN A VIRTUAL ENVIRONMENT FOR EACH USER
            # TODO-ROB The virtual environment can be the name of the user
            # TODO-ROB When the user logs in, they will activate the virtual environment
            for path in self.user_dict.values():
                Path.mkdir(path, parents=True, exist_ok=True)
            for path in self.user_project_dict.values():
                Path.mkdir(path, parents=True, exist_ok=True)

    def user_project_config(self, project_name, private=False):
        if private is False:
            path = self.user_project_dict['public']
        else:
            path = self.user_project_dict['private']
        Path.mkdir(path / Path(project_name))





    def get_newick_dir_map(self, top, ignore=None):
        """Takes a treelib tree created by get_dir_map and returns
        a tree a dir_map in newick format.  This will be useful for Bio.Phylo
        applications."""

        tree = self.get_dir_map(top, ignore)
        Ntree = _parse_json(tree.to_newick_json())
        return Ntree

    def proj_mana(self, js_file):
        with open(js_file, 'r') as file:
            pdr = json.load(file)
        return pdr

    def top_level_struct(self):
        # Create variables for the default directories.
        self.docs = self.__project_home / Path('docs')
        self.lib = self.__project_home / Path('lib')
        self.misc = self.__project_home / Path('misc')
        self.users = self.__project_home / Path('users')
        self.__user_home = self.users / Path(os.environ['USER'])  # TODO-ROB:  Set 'USER' environment variable w/ flask
        # TODO-ROB:  configure during registration for all usernames to be lowercase
        self.web = self.__project_home / Path('web')

    def project(self, p_rt):
        """The projects function helps keep up with different projects.
        Create a different function for every new project."""
        # p_rt is the project research target dictionary
        # the key is the project name and the value is a list of target areas for research
        # Example:  {"GPCR_Project" : ["comparative_evolution", "comparative_polymorphism", "natural_selection"]}
        self.orthologs()

        for project in projects:
            Path.mkdir(self.__project_home / Path(project))


    def orthologs(self):
        """Orthologs-Project"""

        # lib directory
        self.archives = self.lib / Path('archives')
        self.bash = self.lib / Path('bash')
        self.index = self.lib / Path('index')
        self.log = self.lib / Path('log')
        self.perl = self.lib / Path('perl')
        self.python = self.lib / Path('python')
        self.py_classes = self.python / Path('classes')
        self.r = self.lib / Path('r')
        self.scripts = self.lib / Path('scripts')
        self.multiprocessing = self.scripts / Path('multiprocessing')
        self.templates = self.lib / Path('templates')
        # raw_data directory
        self.gi_lists = self.raw_data / Path('gi_lists')

    def gpcr(self):
        """GPCR-Orthologs-Project"""
        self.APP = self.__project_home
        self.WWW = self.APP / Path('www')
        self.APP_DATA = self.APP / Path('Data')
        self.APP_TESTS = self.APP / Path('App-Tests')
        self.ROB_TEST = self.APP_TESTS / Path('Rob')
        self.SHAE_TEST = self.APP_TESTS / Path('Shaurita')
        self.OTHER_TEST = self.APP_TESTS / Path('Other')

    def create_dirs(self):
        for key, value in self.__dict__.items():
            if str(type(value)).__contains__('pathlib') is True:
                Path.mkdir(value, parents=True, exist_ok=True)

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









