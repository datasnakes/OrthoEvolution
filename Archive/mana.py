import os
from pathlib import Path

import Tools
import ete3
from cookiecutter.hooks import run_script
from cookiecutter.main import cookiecutter

# TODO-ROB once this is a pypi package all of these will be unnecessary
from Datasnakes import Cookies, Orthologs
from Datasnakes import Manager
from Datasnakes.Tools.utils import Tree


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
        # TODO-ROB ADD a REPOsitory destination path (an output directory for cookiecutter)
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

    def create_repo(self):
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
        cookiecutter(str(self.repo_cookie), no_input=no_input, extra_context=e_c, output_dir=self.file_home)

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
    def path_list_make(self, path, o_path=None):
        # Takes a path and reduces it to a list of directories within the project
        # An optional attribute (o_path) is give so that a deeper path within the project can be used
        home = str(self.__project_home).split('/')
        path_list = str(path).split('/')
        for item in home:
            if item in path_list:
                path_list.remove(item)
        # path_list = set(p) - set(home)
        if o_path is not None:
            o_path = str(o_path).split('/')
            for item in o_path:
                if item in path_list:
                    path_list.remove(item)
            # path_list = set(path_list) - set(o_path)
        return path_list

    # //TODO-ROB utilize Path.mkdir(parents=TRUE) instead
        # DEPRECATED Change this IN OTHER CLASSES
    def dir_make(self, path, path_list):
        # Takes a path list which is a list of folder names
        # path_list created by str(path).split('/')
        # The path_list appends to path, which is already an established directory inside the project
        t = None
        for item in path_list:

            if os.path.isdir(path + '/' + item): # If for some reason the directory already exists...
                path += '/' + item  # Append a directory
                continue
            path += '/' + item  # Append a directory
            os.mkdir(path)
        if len(os.listdir(path)) > 0:
            path, t = self.dir_archive(path, path_list='')
        return path, t

    # # //TODO-ROB Change to using a compression module https://pymotw.com/2/compression.html
        # DEPRECATED Change this IN OTHER CLASSES
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

# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class RepoMana(Mana):

    def __init__(self, repo, user=None, home=os.getcwd(), new_user=False, new_repo=False):
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

    def create_user(self):
        """This function uses the username given by our FLASK framework
        and creates a new directory system for the active user using
        our  new_user cookiecutter template."""
        # This is used ONLY when the user registers in flask
        # TODO-ROB:  Create the cookiecutter.json file
        # extra_context overrides user and default configs
        cookiecutter(self.user_cookie, no_input=True, extra_context={"user_name": self.user}, output_dir=self.users)
        # TODO-ROB do we need create user hooks?

# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO-ROB:  Edit the setup.py file for cookiecutter.


class UserMana(RepoMana):
    # TODO-ROB CREATE THESE IN A VIRTUAL ENVIRONMENT FOR EACH USER
    # TODO-ROB The virtual environment can be the name of the user
    # TODO-ROB When the user logs in, they will activate the virtual environment
    # TODO-ROB USE SQL here to see if the user db contains the username
    def __init__(self, repo, user, project=None, home=os.getcwd(), new_user=False, new_project=False):
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
        super().__init__(repo=repo, user=user, home=home, new_user=new_user)
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
        cookiecutter(self.project_cookie, extra_context=e_c, no_input=no_input, output_dir=self.project_path)

    def zip_data(self):
        print('zip the users data and send it to their email')

# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class WebMana(RepoMana):

    def __init__(self, repo, website, host='0.0.0.0', port='5252', home=os.getcwd(), new_website=False, create_admin=False):
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
        super().__init__(repo=repo, home=home)
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
        cookiecutter(str(self.website_cookie), no_input=True, extra_context=e_c, output_dir=self.flask)
        # Get the absolute path to the script that starts the flask server
        script_path = self.website_path / Path('hooks') / Path('post_gen_project.sh')
        #scripts_file_path = find_hook('post_gen_project.sh', hooks_dir=str(script_path))
        # TODO-ROB add screening to the bash script for flask run -h -p
        run_script(script_path=str(script_path), cwd=str(self.website_path))


# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# datasnakes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ProjMana(UserMana):

    def __init__(self, repo, user, project, research=None, research_type=None, app=None, home=os.getcwd(),
                 new_project=False, new_research=False, new_app=False):
        """
        :param repo (string):  The name of the repository.
        :param user (string):  The name of the current user if any.
        :param project(string):  The name of the current project if any.
        :param research (string):  The name of the current type of research if any
        :param research_type (string):  The type of research (public or private)
        :param app (string):  The name of the application that the research.
        :param home (string or pathlike):  The home path of the repository.
        :param new_user (bool):  Flag for creating a new user.
        :param new_project (bool):  Flag for creating a new project.
        :param new_research:
        :param new_app:

        """
        super().__init__(repo=repo, user=user, project=project, home=home,
                         new_project=new_project)
        # TODO-ROB Go back to the drawing board for the public/private/other choices.  (FLASK forms)
        # TODO-ROB determine how to get cookiecutter to skip over directories that already exist
        self.project = project
        self.research = research
        self.research_type = research_type

        # Project/Research Directories
        self.research_path = self.project_path / Path(research_type) / Path(research)
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
        '''
        :param new_app (bool):  Flag for auto generating an app that
         goes with the research target.
        :return:  Adds new directories in the current project labeled
        with the proper names.
        '''
        e_c = {"research_type": self.research_type,
               "research_name": self.research}
        cookiecutter(self.research_cookie, no_input=True, extra_context=e_c, output_dir=self.research_path)
        if new_app is True:
            self.create_app()

    def create_app(self):
        e_c = {"app_name": self.app}
        cookiecutter(self.app_cookie, no_input=True, extra_context=e_c, output_dir=self.app_path)


