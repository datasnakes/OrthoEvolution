"""Directory management tools for the package."""
import os
import pkg_resources
from pathlib import Path
from OrthoEvol import Cookies, Orthologs, Manager, Tools
from OrthoEvol.Cookies import Oven
from OrthoEvol.Tools.zipper.zipper import ZipUtils
from OrthoEvol.Tools.logit import LogIt


class Management(object):

    def __init__(self, repo=None, home=os.getcwd(), new_repo=False, **kwargs):
        """
        This is a base class for directory management.

        It maps the directories of the OrthoEvol-Script package using the pathlib module, and turns the names of each
        important directory into a pathlike object.  The base class gives the option of creating a new repository with
        cookiecutter.

        :param repo(string): The name of the new repository to be created.
        :param home(path or path-like): The home of the file calling this name.  When creating a new
            repository it is best to explicitly name the home path.
        :param new_repo(bool): Triggers cookiecutter to create a new repository.
        """

        self.repo = repo
        self.file_home = Path(home)  # Home of the file calling this class
        self.managementlog = LogIt().default(logname="Management", logfile=None)

        # Below are path-like attributes that map various modules and directories.
        # Cookies Module:
        self.Kitchen = Oven(repo=self.repo, output_dir=self.file_home)
        self.Pantry = self.Kitchen.Recipes
        # Manager Module:
        self.Manager = Path(pkg_resources.resource_filename(Manager.__name__, ''))
        self.BioSQL = self.Manager / Path('BioSQL')
        self.SQLite3 = self.BioSQL / Path('sqlite')
        self.MySQL = self.BioSQL / Path('mysql')
        self.config = self.Manager / Path('config')
        # Orthologs Module:
        self.Orthologs = Path(pkg_resources.resource_filename(Orthologs.__name__, ''))
        self.Align = self.Orthologs / Path('Align')
        self.Blast = self.Orthologs / Path('Blast')
        self.GenBank = self.Orthologs / Path('GenBank')
        self.Phylogenetics = self.Orthologs / Path('Phylogenetics')
        # Tools Module:
        self.Tools = Path(pkg_resources.resource_filename(Tools.__name__, ''))
        self.ftp = self.Tools / Path('ftp')
        self.logit = self.Tools / Path('logit')
        self.mpi = self.Tools / Path('mpi')
        self.mygene = self.Tools / Path('mygene')
        self.pandoc = self.Tools / Path('pandoc')
        self.parallel = self.Tools / Path('parallel')
        self.pybasher = self.Tools / Path('pybasher')
        self.send2server = self.Tools / Path('send2server')
        self.sge = self.Tools / Path('sge')
        self.slackify = self.Tools / Path('slackify')
        self.utils = self.Tools / Path('utils')
        self.zipper = self.Tools / Path('zipper')

        if self.repo:
            self.repo_path = self.file_home / Path(self.repo)
        self.managementlog.info('The BaseManagement class variables have been set.')

        # Make a new repository.
        if new_repo is True:
            self.managementlog.info('The repository cookie is being prepared for the Oven.')
            self.Kitchen.bake_the_repo()


class RepoManagement(Management):

    def __init__(self, repo, user=None, home=os.getcwd(),
                 new_user=False, new_repo=False, **kwargs):
        """
        This is the Repository Management class, which inherits the Management base class.  This class has to be paired
        with a repository name.  It gives the option of creating a filesystem for a new user and a new repository
        using cookiecutter.

        This class maps the named repository, which includes top level access to front facing web servers, important
        documents, and the top level users directory.

        :param repo (string):  The name of the repository.
        :param user (string):  The name of the current user if any.
        :param home (string or pathlike):  The home path of the repository.
        :param new_user (bool):  Flag for creating a new user.
        :param new_repo (bool): Flag for creating a new repository.
        """

        # TODO-ROB change the home parameter to the output directory parameter
        super().__init__(repo=repo, home=home, new_repo=new_repo, **kwargs)

        self.repo = repo
        # Users and Important Documentation:
        self.docs = self.repo_path / Path('docs')
        self.users = self.repo_path / Path('users')
        # Web Servers:
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
            self.managementlog.info('The user cookie is being prepared for the Oven.')
            self.Kitchen.bake_the_user()
# TODO-ROB:  Edit the setup.py file for cookiecutter.


class UserManagement(RepoManagement):
    # TODO-ROB CREATE THESE IN A VIRTUAL ENVIRONMENT FOR EACH USER
    # TODO-ROB The virtual environment can be the name of the user
    # TODO-ROB When the user logs in, they will activate the virtual environment
    # TODO-ROB USE SQL here to see if the user db contains the username

    def __init__(self, repo, user, project=None, db_config_file=None, home=os.getcwd(),
                 new_user=False, new_project=False, new_db=False, archive=False, **kwargs):
        """
        This is the User Management class, which manages the current users directories.  This class has to be paired
        with a repository and a user.  It gives access to user paths, and provides functionality for creating new
        projects and new project db_config_file for the current user.  It also give the option of creating a new user.

        This class maps a users directory, which gives access to directories for db_config_file (NCBI and proprietary), index
        files for quickly retrieving project data, project log files, user affiliated journal articles, and projects.

        :param repo (string):  The name of the repository.
        :param user (string):  The name of the current user if any.
        :param project(string):  The name of the current project if any.
        :param home (string or pathlike):  The home path of the repository.
        :param new_user (bool):  Flag for creating a new user.
        :param new_project (bool):  Flag for creating a new project.
        """

        if user or (user and repo):
            super().__init__(repo=repo, user=user, home=home, new_user=new_user, **kwargs)
            self.user = user
            # NCBI and Proprietary Database Repositories:
            self.user_db = self.user_path / Path('databases')
            self.user_archive = self.user_path / Path('archive')
            self.ncbi_db_repo = self.user_db / Path('NCBI')
            self.ncbi_taxonomy = self.ncbi_db_repo / Path('pub') / Path('taxonomy')
            self.ncbi_refseq_release = self.ncbi_db_repo / Path('refseq') / Path('release')
            self.blast_db = self.ncbi_db_repo / Path('blast') / Path('db')
            self.windowmaker_files = self.ncbi_db_repo / Path('blast') / Path('windowmaker_files')
            self.itis_db_repo = self.user_db / Path('ITIS')
            # Index Files:
            self.user_index = self.user_path / Path('index')
            # User Log Files:
            self.user_log = self.user_path / Path('log')
            # Relevant Journal Articles:
            self.manuscripts = self.user_path / Path('manuscripts')
            self.other = self.user_path / Path('other')
            # Projects
            self.projects = self.user_path / Path('projects')

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

        self.managementlog.info('The User Management class variables have been set.')

        if new_project is True:
            self.managementlog.info('The project cookie is being prepared for the Oven.')
            self.Kitchen.bake_the_project()
        if new_db is True:
            self.managementlog.info('The database cookie is being prepared for the Oven.')
            self.Kitchen.bake_the_db_repo(user_db=self.user_db, db_path_dict=self.db_path_dict, ncbi_db_repo=self.ncbi_db_repo)
            # TODO-ROB:  Determine what type of database as well.

    def zip_mail(self, comp_filename, zip_path, destination=''):
        Zipper = ZipUtils(comp_filename, zip_path)
        Zipper_path = Zipper.compress()
        # TODO-ROB add proper destination syntax.
        self.managementlog.info('%s is being sent to %s' % (Zipper_path, destination))


class WebsiteManagement(RepoManagement):

    def __init__(self, repo, website, host='0.0.0.0', port='5252',
                 home=os.getcwd(), new_website=False, create_admin=False, **kwargs):
        """This is the Website Management class, which installs a template for Flask using cookiecutter.  The
        official cookiecutter-flask template (https://github.com/sloria/cookiecutter-flask) has been edited for our own
        purposes.  This app class uses cookiecutter hooks to deploy the flask server.

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

        # Website Deployment Information:
        self.website = website
        self.web_host = host
        self.web_port = port
        # Path to Flask's Web-Server Files
        self.website_path = self.flask / Path(self.website)

        self.Kitchen = Oven(repo=self.repo, user=self.user, website=self.website, output_dir=self.flask)
        self.managementlog.info('The Website Management class variables have been set.')

        if new_website is True:
            self.managementlog.info('The website cookie is being prepared for the Oven.')
            self.Kitchen.bake_the_website(host=self.web_host, port=self.web_port, website_path=self.website_path)

    def stop_server(self):
        """Stop the server running the website."""
        # TODO-SDH Add way to stop the server from running.


class ProjectManagement(UserManagement):

    def __init__(self, repo, user, project, research=None, research_type=None,
                 app=None, home=os.getcwd(), new_project=False, new_research=False,
                 new_app=False, **kwargs):
        """
        This is the Project Management class, which manages the directories of the current project.  Each project
        requires a repository, user, and project name.  It gives the option of starting a new type of research within
        an existing project.  An application directory for the specific research/dataset can also be generated

        It gives access to the project directories including index
        files, the raw data, the processed data, the project db_config_file, and the web files for serving data.

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

            self.project = project
            self.research = research
            self.research_type = research_type
            # Project Directories:
            self.research_path = self.project_path / Path(research_type) / Path(research)
            self.project_archive = self.project_path / Path('archive')
            self.project_database = self.user_db / Path(project)
            # Dataset Directories:
            self.project_index = self.research_path / Path('index')
            self.data = self.research_path / Path('data')
            self.raw_data = self.research_path / Path('raw_data')
            self.project_web = self.research_path / Path('web')
            if app:
                self.app = app
                self.app_path = self.project_web / Path(app)

        self.managementlog.info('The User Management class variables have been set.')

        if new_research is True:
            self.managementlog.info('The research cookie is being prepared for the Oven.')
            self.research_type = research_type
            self.Kitchen = Oven(repo=self.repo, user=self.user, project=self.project, output_dir=self.project_path)
            self.Kitchen.bake_the_research(research_type=self.research_type, research=self.research)
            if new_app is True:
                self.managementlog.info('The app cookie is being prepared for the Oven.')
                self.app = app
                self.app_path = self.project_path / Path(research_type) / Path(research) / Path('web')
                self.Kitchen.bake_the_app(app=self.app)

