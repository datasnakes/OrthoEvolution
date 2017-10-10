"""Management tools for the package."""
import os
import pkg_resources
from pathlib import Path
from Datasnakes import Cookies, Orthologs, Manager, Tools
from Datasnakes.Cookies.cookie_jar import Oven
from Datasnakes.Tools.zipper.zipper import ZipUtils
from Datasnakes.Tools.logit import LogIt


class Management(object):
    """This is the directory management base class.

    It maps the directories in the PyPi package using the pathlib module and
    turns the names of each important directory into a pathlike object.  The
    base class gives the option of creating a new repository with cookiecutter.

    This is also the home for many of the utility functions for manipulating
    directories or paths.
    """

    def __init__(self, repo=None, home=os.getcwd(), new_repo=False, **kwargs):
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
        self.managementlog = LogIt().default(logname="Management", logfile=None)
        # TODO-ROB:  SOme of these directories don't need to be accessed directly
        # Below are the PyPi path strings
        #    The first group is to access the cookiecutter templates
        self.Kitchen = Oven(repo=self.repo, output_dir=self.file_home)
        self.Pantry = self.Kitchen.Ingredients
        #    The second group is for the Manager module
        self.Manager = Path(pkg_resources.resource_filename(Manager.__name__, ''))
        self.BioSQL = self.Manager / Path('BioSQL')
        self.config = self.Manager / Path('config')
        #    The third group is for the Orthologs module
        self.Orthologs = Path(pkg_resources.resource_filename(Orthologs.__name__, ''))
        self.Align = self.Orthologs / Path('Align')
        self.Blast = self.Orthologs / Path('Blast')
        self.GenBank = self.Orthologs / Path('GenBank')
        self.Phylogenetics = self.Orthologs / Path('Phylogenetics')
        #    The fourth group is for the Tools module
        self.Tools = Path(pkg_resources.resource_filename(Tools.__name__, ''))
        self.ftp = self.Tools / Path('ftp')
        self.logit = self.Tools / Path('logit')
        self.mpi = self.Tools / Path('mpi')
        self.mygene = self.Tools / Path('mygene')
        self.pandoc = self.Tools / Path('pandoc')
        self.parallel = self.Tools / Path('parallel')
        self.pybasher = self.Tools / Path('pybasher')
        self.qsub = self.Tools / Path('qsub')
        self.send2server = self.Tools / Path('send2server')
        self.slackify = self.Tools / Path('slackify')
        self.utils = self.Tools / Path('utils')
        self.zipper = self.Tools / Path('zipper')

        if self.repo:
            self.repo_path = self.file_home / Path(self.repo)

        self.managementlog.info('The BaseManagement class variables have been set.')

        if new_repo is True:
            self.Kitchen.bake_the_repo()



        # Create a directory management logger
        # TODO-ROB add logging to manager class
        #log = LogIt('user/path/userfile.log', 'Directory Management')
        #self.dm_log = log.basic


class RepoManagement(Management):
    """Repository Management."""
    def __init__(self, repo, user=None, home=os.getcwd(),
                 new_user=False, new_repo=False, **kwargs):
        """
        :param repo (string):  The name of the repository.
        :param user (string):  The name of the current user if any.
        :param home (string or pathlike):  The home path of the repository.
        :param new_user (bool):  Flag for creating a new user.
        :param new_repo: Flag for creating a new repository.
        """
        # TODO-ROB change the home parameter to the output directory parameter
        super().__init__(repo=repo, home=home, new_repo=new_repo, **kwargs)
        self.repo = repo
        self.docs = self.repo_path / Path('docs')
        self.users = self.repo_path / Path('users')

        self.repo_web = self.repo_path / Path('web')
        self.repo_shiny = self.repo_web / Path('shiny')
        self.ftp = self.repo_web / Path('ftp')
        self.wasabi = self.repo_web / Path('wasabi')
        self.flask = self.repo_web / Path('flask')

        if user:
            self.user = user  # FROM Flask
            self.user_path = self.users / Path(self.user)

        self.Kitchen = Oven(repo=self.repo, user=self.user, output_dir=self.users)
        self.managementlog.info('The Repository Management class variables have been set.')

        if new_user is True:
            self.Kitchen.bake_the_user()
        # TODO-ROB do we need create user hooks?
# TODO-ROB:  Edit the setup.py file for cookiecutter.


class UserManagement(RepoManagement):
    """User Management Class."""
    # TODO-ROB CREATE THESE IN A VIRTUAL ENVIRONMENT FOR EACH USER
    # TODO-ROB The virtual environment can be the name of the user
    # TODO-ROB When the user logs in, they will activate the virtual environment
    # TODO-ROB USE SQL here to see if the user db contains the username

    def __init__(self, repo, user, project=None, database=None, home=os.getcwd(),
                 new_user=False, new_project=False, new_db=False, **kwargs):
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
        if database is None:
            database = []
        if user or (user and repo):
            super().__init__(repo=repo, user=user, home=home, new_user=new_user, **kwargs)
            self.user = user

            # TODO-ROB Add database files to the repository as well
            self.user_db = self.user_path / Path('databases')
            self.ncbi_db_repo = self.user_db / Path('NCBI')
            self.user_index = self.user_path / Path('index')
            self.user_log = self.user_path / Path('log')
            self.manuscripts = self.user_path / Path('manuscripts')
            self.other = self.user_path / Path('other')
            self.projects = self.user_path / Path('projects')
        # TODO-ROB Create a DB_mana class in a seperate file that interactsd with the ftp class
            if project:
                self.project = project
                self.project_path = self.projects / Path(project)
        else:
            self.projects = home
            self.Cookies = Path(Cookies.__path__[0])
            self.project_cookie = self.Cookies / Path('new_project')
            if project:
                self.project = project
                self.project_path = home / Path(project)
        self.Kitchen = Oven(repo=self.repo, user=self.user, project=self.project, output_dir=self.projects)

        if len(database) > 0:
            self.db_list = database
            self.db_path_dict = {}
            self.db_archive_dict = {}
            for item in database:
                self.db_path_dict[item] = self.user_db / Path(item)
                self.db_archive_dict[item] = self.db_path_dict[item] / Path('archive')
        else:
            self.db_path_dict = None
            self.db_archive_dict = None

        self.managementlog.info('The User Management class variables have been set.')

        if new_project is True:
            self.Kitchen.bake_the_project()
        if new_db is True:
            self.Kitchen.bake_the_db_repo(user_db=self.user_db, db_path_dict=self.db_path_dict, ncbi_db_repo=self.ncbi_db_repo)

    def zip_mail(self, comp_filename, zip_path, destination=''):
        Zipper = ZipUtils(comp_filename, zip_path)
        Zipper_path = Zipper.to_zip()
        # TODO-ROB add proper destination syntax.
        self.managementlog.info('%s is being sent to %s' % (Zipper_path, destination))


class WebsiteManagement(RepoManagement):
    """Web Management Class."""

    def __init__(self, repo, website, host='0.0.0.0', port='5252',
                 home=os.getcwd(), new_website=False, create_admin=False, **kwargs):
        """Install a template for Flask using cookiecutter.

        The custom datasnakes cookie for this template has been edited for
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
        """
        if repo and website:
            super().__init__(repo=repo, home=home, **kwargs)
        self.website = website
        self.web_host = host
        self.web_port = port
        self.website_path = self.flask / Path(self.website)

        self.Kitchen = Oven(repo=self.repo, user=self.user, website=self.website, output_dir=self.flask)

        self.managementlog.info('The Website Management class variables have been set.')

        if new_website is True:
            self.Kitchen.bake_the_website(host=self.web_host, port=self.web_port, website_path=self.website_path)

    def stop_server(self):
        """Stop the server running the website."""
        # TODO-SDH Add way to stop the server from running.


class ProjectManagement(UserManagement):
    """Project Management Class."""

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
        # Standalone for child/self or full class hierarchy use
        if project or (repo and user and project):
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
            self.project_database = self.user_db / Path(project)
            self.data = self.research_path / Path('data')
            self.raw_data = self.research_path / Path('raw_data')
            self.project_web = self.research_path / Path('web')
            # TODO-ROB:  THis is just a draft.  Rework to use public/private/other
            if app:
                self.app = app
                self.app_path = self.project_web / Path(app)

        self.managementlog.info('The User Management class variables have been set.')

        if new_research is True:
            self.research_type = research_type
            self.Kitchen = Oven(repo=self.repo, user=self.user, project=self.project, output_dir=self.project_path)
            self.Kitchen.bake_the_research(research_type=self.research_type, research=self.research)
            if new_app is True:
                self.app = app
                self.app_path = self.project_path / Path(research_type) / Path(research) / Path('web')
                self.Kitchen.bake_the_app(app=self.app)

