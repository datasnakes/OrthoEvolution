# -*- coding: utf-8 -*-
# Modules Used
import os
import stat
from pathlib import Path
from Manager.utils.zipper import ZipUtils
from cookiecutter.main import cookiecutter
from cookiecutter.hooks import run_script, find_hook
from Manager.utils.treelib2.treelib2.tree import Tree
import ete3
from Manager.utils.json_to_newick import _parse_json
# TODO-ROB once this is a pypi package all of these will be unnecessary
import Cookies
import Manager
import Orthologs
import Tools
# from Manager.logit.logit import LogIt

#------------------------------------------------------------------------------
# TODO-ROB use **kwargs and **args to cut down on parameters


class Mana(object):
    """
    This is the directory management base class.  It
    maps the directories in the PyPi package using the pathlib module and
    turns the names of each important directory into a pathlike object.  The
    base class gives the option of creating a new repository with cookiecutter.

    This is also the home for many of the utility functions for manipulating
    directories or paths.
    """

    def __init__(self, repo=None, home=os.getcwd(), new_repo=False):
        # TODO-ROB ADD a REPOsitory destination path (an output directory for
        # cookiecutter)
        """
        :param home(path or path-like): The home of the file calling this name.  When creating a new
            repository it is best to explicitly name the home path.
        :param repo(string): The name of the new repository to be created.
        :param new_repo(bool): Triggers cookiecutter to create a new repository.
        """
        self.repo = repo
        self.file_home = Path(home)  # Home of the file calling this class
        # TODO-ROB:  SOme of these directories don't need to be accessed directly
        # Below are the PyPi path strings
        #    The first group is to access the cookiecutter templates
        self.Cookies = Path(Cookies.__path__._path[0])
        self.repo_cookie = self.Cookies / Path('new_repository')
        self.user_cookie = self.Cookies / Path('new_user')
        self.project_cookie = self.Cookies / Path('new_project')
        self.research_cookie = self.Cookies / Path('new_research')
        self.app_cookie = self.Cookies / Path('new_app')
        self.website_cookie = self.Cookies / Path('new_website')
        #    The second group is for the Manager module
        self.Manager = Path(Manager.__path__._path[0])
        self.index = self.Manager / Path('index')
        self.logit = self.Manager / Path('logit')
        self.utils = self.Manager / Path('utils')
        self.shiny = self.Manager / Path('shiny')
        #    The third group is for the Orthologs module
        self.Orthologs = Path(Orthologs.__path__._path[0])
        self.biosql = Path(self.Orthologs) / Path('biosql')
        self.blast = Path(self.Orthologs) / Path('blast')
        self.comp_gen = Path(self.Orthologs) / Path('comparative_genetics')
        self.genbank = Path(self.Orthologs) / Path('genbank')
        self.manager = Path(self.Orthologs) / Path('manager')
        self.phylogenetics = Path(self.Orthologs) / Path('phylogenetics')
        #    The fourth group is for the Tools module
        self.Tools = Path(Tools.__path__._path[0])
        self.ftp = Path(self.Tools) / Path('ftp')
        self.multiprocessing = Path(self.Tools) / Path('multiprocessing')
        self.pandoc = Path(self.Tools) / Path('pandoc')
        self.pybasher = Path(self.Tools) / Path('pybasher')
        self.qsub = Path(self.Tools) / Path('qsub')

        if self.repo:
            self.repo_path = self.file_home / Path(self.repo)
        if new_repo is True:
            self.create_repo()

        # Create a directory management logger
        # TODO-ROB add logging to manager class
        #log = LogIt('user/path/userfile.log', 'Directory Management')
        #self.dm_log = log.basic

#------------------------------------------------------------------------------
    def create_repo(self):
        print('creating dirs from repo cookie')
        print(self.__class__.__name__)
        """This function creates a new repository.  If a repository name
        is given to the class then it is given a name.  If not, cookiecutters
        takes input from the user.

        The base class will be the only class that allows cookiecutters parameter
        no_input to be False.
        """
        if self.repo:
            no_input = True
            e_c = {
                "repository_name": self.repo
            }
        else:
            no_input = False
            e_c = None
            # TODO-ROB change cookiecutter so that it can take pathlike objects
        cookiecutter(str(self.repo_cookie), no_input=no_input,
                     extra_context=e_c, output_dir=str(self.file_home))
        os.chmod(str(self.file_home / Path(self.repo)), mode=0o777)

#------------------------------------------------------------------------------
    # def git_ignore(self, path):
    #     """Get the ignored file patterns from the .gitignore file in the repo."""
    #     with open(Path(path) / Path('.gitignore'), 'r', newline='') as ignore:
    #         ignored = ignore.read().splitlines()
    #     return ignored
    #
    # # Map the main project directory.
    # def get_dir_map(self, top, gitignore=None):
    #     # TODO-ROB:  Change ignore to a .gitignore filename
    #     default_ignore = self.git_ignore(top)
    #     if gitignore is not None:
    #         gitignore += default_ignore
    #     else:
    #         gitignore = default_ignore
    #     # Treelib will help to map everything and create a json at the same time
    #     tree = Tree()
    #     tree.create_node('.', top)
    #     # Walk through the directory of choice (top)
    #     # Topdown is true so that directories can be modified in place
    #     for root, dirs, files in os.walk(top, topdown=True):
    #         # Only remove directories from the top
    #         if root == top:
    #             print(root)
    #             try:
    #                 dirs[:] = set(dirs) - set(gitignore)  # Remove directories from os.walk()
    #                 print(dirs)
    #             except ValueError:
    #                 pass
    #         for d in dirs:
    #             rd = str(Path(root) / Path(d))
    #             tree.create_node(d, identifier=rd, parent=root)
    #         for f in files:
    #             tree.create_node(f, parent=root)
    #     return tree

#------------------------------------------------------------------------------
    def get_newick_dir_map(self, top, ignore=None):
        """Takes a treelib tree created by get_dir_map and returns
        a tree a dir_map in newick format.  This will be useful for Bio.Phylo
        applications."""
        """
        :param top (path):  The root at which the directory map is made.
        :param ignore (list):  The files to ignore.  The  get_dir_map function
        adds this to the .gitignore list.
        :return (tree):  A newick formatted string in style #8.  Can be used with
        the ete3.Tree() class.
        """

        tree = Tree()
        t = tree.get_dir_map(top, ignore)
        Ntree = tree.parse_newick_json()
        return Ntree

    def get_ete3_tree(self, top, tree=None):
        if not tree:
            tree = self.get_newick_dir_map(top)
        t = ete3.Tree(tree, format=8)
        return t

    # DEPRECATED Change this IN OTHER CLASSES
   # def path_list_make(self, path, o_path=None):
        # Takes a path and reduces it to a list of directories within the project
        # An optional attribute (o_path) is give so that a deeper path within
        # the project can be used

    # //TODO-ROB utilize Path.mkdir(parents=TRUE) instead
        # DEPRECATED Change this IN OTHER CLASSES
  #  def dir_make(self, path, path_list):
        # Takes a path list which is a list of folder names
        # path_list created by str(path).split('/')
        # The path_list appends to path, which is already an established
        # directory inside the project

    # # //TODO-ROB Change to using a compression module https://pymotw.com/2/compression.html
        # DEPRECATED Change this IN OTHER CLASSES
    # def dir_archive(self, path, path_list):
    #     # Use the path that you want to update/add to
    #     # Returns path and the time stamp (could be None)

#------------------------------------------------------------------------------


class RepoMana(Mana):

    def __init__(self, repo, user=None, home=os.getcwd(),
                 new_user=False, new_repo=False):
        """
        :param repo (string):  The name of the repository.
        :param user (string):  The name of the current user if any.
        :param home (string or pathlike):  The home path of the repository.
        :param new_user (bool):  Flag for creating a new user.
        :param new_repo: Flag for creating a new repository.
        """
        # TODO-ROB change the home parameter to the output directory parameter
        super().__init__(repo=repo, home=home, new_repo=new_repo)
        self.repo = repo
        self.docs = self.repo_path / Path('docs')
        self.misc = self.repo_path / Path('misc')
        self.users = self.repo_path / Path('users')

        self.lib = self.repo_path / Path('lib')
        self.archives = self.lib / Path('archives')
        self.databases = self.lib / Path('databases')
        self.repo_index = self.lib / Path('index')

        self.repo_web = self.repo_path / Path('web')
        self.repo_shiny = self.repo_web / Path('shiny')
        self.ftp = self.repo_web / Path('ftp')
        self.wasabi = self.repo_web / Path('wasabi')

        self.flask = self.repo_web / Path('flask')

        if user:
            self.user = user  # FROM Flask
            self.user_path = self.users / Path(self.user)
        if new_user is True:
            self.create_user()

#------------------------------------------------------------------------------
    def create_user(self):
        """This function uses the username given by our FLASK framework
        and creates a new directory system for the active user using
        our  new_user cookiecutter template.
        """
        print('creating dirs from user cookie')
        print(self.__class__.__name__)

        # This is used ONLY when the user registers in flask
        # TODO-ROB:  Create the cookiecutter.json file

        # extra_context overrides user and default configs
        cookiecutter(str(self.user_cookie), no_input=True, extra_context={
                     "user_name": self.user}, output_dir=str(self.users))

        # Change user permissions with flask later (this is for testing
        # purposes
        os.chmod(str(self.users / Path(self.user)), mode=0o777)
        # TODO-ROB do we need create user hooks?

# TODO-ROB:  Edit the setup.py file for cookiecutter.
#------------------------------------------------------------------------------


class UserMana(RepoMana):
    """User Management Class.
    """
    # TODO-ROB CREATE THESE IN A VIRTUAL ENVIRONMENT FOR EACH USER
    # TODO-ROB The virtual environment can be the name of the user
    # TODO-ROB When the user logs in, they will activate the virtual environment
    # TODO-ROB USE SQL here to see if the user db contains the username

    def __init__(self, repo, user, project=None, home=os.getcwd(),
                 new_user=False, new_project=False, **kwargs):
        '''
        The User Management class manages the current users directories.
        This class gives access to user paths, and provides functionality
        for creating new projects for the current user within the users
        home.

        :param repo (string):  The name of the repository.
        :param user (string):  The name of the current user if any.
        :param project(string):  The name of the current project if any.
        :param home (string or pathlike):  The home path of the repository.
        :param new_user (bool):  Flag for creating a new user.
        :param new_project (bool):  Flag for creating a new project.
        '''
        super().__init__(repo=repo, user=user, home=home, new_user=new_user, **kwargs)
        self.user = user

        self.user_index = self.user_path / Path('index')
        self.user_log = self.user_path / Path('log')
        self.manuscripts = self.user_path / Path('manuscripts')
        self.other = self.user_path / Path('other')
        self.projects = self.user_path / Path('projects')

        if project:
            self.project = project
            self.project_path = self.projects / Path(project)
        if new_project is True:
            self.create_project()

    def create_project(self):
        print('creating dirs from project cookie')
        print(self.__class__.__name__)
        """
        :return: A new project inside the user's
        project directory.
        """
        if self.project:
            no_input = True
            e_c = {"project_name": self.project}
        else:
            no_input = False
            e_c = None
        cookiecutter(str(self.project_cookie), extra_context=e_c,
                     no_input=no_input, output_dir=str(self.projects))
        os.chmod(str(self.projects / Path(self.project)), mode=0o777)

    def zip_mail(self, comp_filename, zip_path, destination=''):
        Zipper = ZipUtils(comp_filename, zip_path)
        Zipper_path = Zipper.to_zip()
        # TODO-ROB add proper destination syntax.
        print('%s is being sent to %s' % (Zipper_path, destination))

#------------------------------------------------------------------------------


class WebMana(RepoMana):
    """Web Management Class.
    """

    def __init__(self, repo, website, host='0.0.0.0', port='5252',
                 home=os.getcwd(), new_website=False, create_admin=False, **kwargs):
        '''
        This installs a template for Flask using cookiecutter.  The
        custom datasnakes cookie for this template has been edited for
        our own purposes.

        :param repo (string):  The name of the repository.
        :param website (string):  The name of the website.  Not a url, so
        it doesn't containt http://www.*.com.
        (e.g. for www.vallenger-genetics.ml this parameter would be 'vallender-genetics')
        :param host (string):  The address to launch the flask app.  Defaults to 0.0.0.0
        :param port (string):  The port to launch the flask app.  Defaults to 5252
        :param home (string or pathlike):  The home path of the repository.
        :param new_website (bool):  Flag for creating a new website
        :param create_admin:  Flag for creating a new admin for the website via FLASK USER.
        (Note:  This parameter is not used currently in development.)
        '''
        super().__init__(repo=repo, home=home, **kwargs)
        self.website = website
        self.web_host = host
        self.web_port = port
        self.website_path = self.flask / Path(self.website)

        self.website_scripts = self.website_path / Path(self.website)
        self.website_public = self.website_scripts / Path('public')
        self.website_user = self.website_scripts / Path('user')

        if new_website is True:
            self.create_website()

    def create_website(self):
        '''
        To create a website, the new_website cookie is used.
        After creating the directory structure, the run_script function
        from cookiecutter finds the hooks folder which contains a
        post-cookiecutter-template-generation bash script.  The bash script
        sets up the proper dependencies and environment variables for the website,
        and runs the website on the specified host and port

        :return: Runs the website.
        '''
        # TODO-ROB Add heavy logging here
        e_c = {"website_name": self.website,
               "website_path": os.path.join(str(self.website_path), ''),
               "website_host": self.web_host,
               "website_port": self.web_port}
        cookiecutter(str(self.website_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.flask))
        os.chmod(str(self.flask / Path(self.website)), mode=0o777)
        # Get the absolute path to the script that starts the flask server
        script_path = self.website_path / \
            Path('hooks') / Path('post_gen_project.sh')
        #scripts_file_path = find_hook('post_gen_project.sh', hooks_dir=str(script_path))
        # TODO-ROB add screening to the bash script for flask run -h -p
        run_script(script_path=str(script_path), cwd=str(self.website_path))

#------------------------------------------------------------------------------


class ProjMana(UserMana):
    """Project Management Class.
    """

    def __init__(self, repo, user, project, research=None, research_type=None,
                 app=None, home=os.getcwd(), new_project=False, new_research=False,
                 new_app=False, **kwargs):
        """
        :param repo (string):  The name of the repository.
        :param user (string):  The name of the current user if any.
        :param project(string):  The name of the current project if any.
        :param research (string):  The name of the current type of research if any
        :param research_type (string):  The type of research (public or private)
        :param app (string):  The name of the application that the research.
        :param home (string or pathlike):  The home path of the repository.
        :param new_project (bool):  Flag for creating a new project.
        :param new_research (bool):  Flag for creating new research under a project.
        :param new_app (bool):  Flag for creating a new web app under a research target.

        """
        super().__init__(repo=repo, user=user, project=project, home=home,
                         new_project=new_project, **kwargs)
        # TODO-ROB Go back to the drawing board for the public/private/other choices.  (FLASK forms)
        # TODO-ROB determine how to get cookiecutter to skip over directories
        # that already exist
        self.project = project
        self.research = research
        self.research_type = research_type
        # Project/Research Directories
        self.research_path = self.project_path / \
            Path(research_type) / Path(research)
        self.project_archive = self.project_path / Path('archive')
        self.project_index = self.research_path / Path('index')
        self.data = self.research_path / Path('data')
        self.raw_data = self.research_path / Path('raw_data')
        self.project_web = self.research_path / Path('web')
        # TODO-ROB:  THis is just a draft.  Rework to use public/private/other
        if app:
            self.app = app
            self.app_path = self.project_web / Path(app)
        if new_research is True:
            self.create_research(new_app)

    def create_research(self, new_app=False):
        """
        :param new_app (bool):  Flag for auto generating an app that
         goes with the research target.
        :return:  Adds new directories in the current project labeled
        with the proper names.
        """
        print('Creating directories from research cookie...')
        print(self.__class__.__name__)

        e_c = {"research_type": self.research_type,
               "research_name": self.research}
        cookiecutter(str(self.research_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.project_path))
        os.chmod(str(self.project_path / Path(self.research_type)), mode=0o777)
        if new_app is True:
            self.create_app()

    def create_app(self):
        e_c = {"app_name": self.app}
        cookiecutter(str(self.app_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.app_path))
        os.chmod(str(self.app_path / Path(self.app)), mode=0o777)
