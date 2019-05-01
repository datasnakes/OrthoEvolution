# Standard Library
import os
from pathlib import Path

# from OrthoEvol.Cookies.utils import archive
# Other
import yaml
from cookiecutter.hooks import run_script
from cookiecutter.main import cookiecutter
from pkg_resources import resource_filename

# OrthoEvol
from OrthoEvol import Cookies
from OrthoEvol.Manager.config import yml
from OrthoEvol.Tools.logit import LogIt


class CookBook(object):
    _config_file = resource_filename(yml.__name__, "cookie_recipes.yml")

    def __init__(self, config_file=_config_file, **new_recipes):
        """
        The Cookie Recipes are public attributes for accessing
        the paths to the various cookiecutter templates in the
        Cookies module.

        The Cookie Recipes are used for the Recipes attribute
        in the Oven class.

        New Recipes can also be added...
        """
        self.CookieJar = Path(resource_filename(Cookies.__name__, ''))
        self.repo_cookie = self.CookieJar / Path('new_repository')
        self.user_cookie = self.CookieJar / Path('new_user')
        self.project_cookie = self.CookieJar / Path('new_project')
        self.basic_project_cookie = self.CookieJar / Path('new_basic_project')
        self.research_cookie = self.CookieJar / Path('new_research')
        self.app_cookie = self.CookieJar / Path('new_app')
        self.db_cookie = self.CookieJar / Path('new_database')
        self.website_cookie = self.CookieJar / Path('new_website')

        # TODO-ROB:  Make this function better.
        # Load the cookies that are in the cookie_jar config file
        with open(config_file, 'r') as ymlfile:
            configuration = yaml.safe_load(ymlfile)
            if configuration is not None:
                setattr(self, "CONFIGURATION", configuration)
                for key, value in configuration.items():
                    setattr(self, key, value)

        # For custom cookies use a dictionary to create attributes
        if new_recipes:
            for cookie, path in new_recipes.items():
                setattr(self, cookie, path)
            # Add the new recipes to the configuration file dictionary
            configuration.update(new_recipes)

        # Overwrite the yaml config file.
        with open(config_file, 'w') as ymlfile:
            yaml.dump(configuration, ymlfile)


class Oven(object):

    def __init__(self, repo=None, user=None, project=None, basic_project=False, website=None, db_repo="databases",
                 output_dir=os.getcwd(), recipes=CookBook()):
        """
        This class uses cookiecutter to deploy custom cookiecutter templates:

        The Oven uses the different Ingredients (parameters/attributes) and
        the Cook Book(cookiecutter templates) to bake_the_cookies
        in the Oven(class methods).

        After the cookies cool, they are put in the cookie_jar (output directory).

        :param repo (string):  An ingredient representing the repository name.
        :param user (string):  An ingredient representing the user name
        :param project (string):  An ingredient representing the project name.
        :param basic_project (bool):  A secret ingredient ONLY for the basic project cookie.
        :param db_config_file (list):  An ingredient representing a list of db_config_file.
        :param website (string):  An ingredient representing the website name.
        :param output_dir (path or pathlike):  The cookie jar for storing the cookies.
        :param recipes (pathlike):  An index for the different recipe templates.
        """
        self.cookielog = LogIt().default(logname="Cookies", logfile=None)
        self.cookie_jar = output_dir
        # Below are the PyPi path strings
        #    The first group is to access the cookiecutter templates
        self.repo = repo
        self.user = user
        self.project = project
        self.basic_project = basic_project
        self.website = website
        self.db_repo = db_repo
        self.Recipes = recipes
        self.Ingredients = {"repo": self.repo,
                            "user": self.user,
                            "project": self.project,
                            "basic_project": self.basic_project,
                            "website": self.website,
                            "db_repo": self.db_repo,
                            "recipes": self.Recipes.__dict__}

    def bake_the_repo(self, cookie_jar=None):
            self.cookielog.warn('Creating directories from the Repository Cookie template.')
            """
            This function creates a new repository.  If a repository name
            is given to the class then it is given a name.  If not, cookiecutters
            takes input from the user.
    
            The base class will be the only class that allows cookiecutters parameter
            no_input to be False.
            """
            if cookie_jar:
                self.cookie_jar = cookie_jar
            if self.repo:
                no_input = True
                e_c = {
                    "repository_name": self.repo
                }
            else:
                no_input = False
                e_c = None
                # TODO-ROB change cookiecutter so that it can take pathlike objects
            cookiecutter(str(self.Recipes.repo_cookie), no_input=no_input,
                         extra_context=e_c, output_dir=str(self.cookie_jar))
            os.chmod(str(self.cookie_jar / Path(self.repo)), mode=0o777)
            self.cookielog.info('Repository directories have been created. ✔')

    def bake_the_user(self, cookie_jar=None):
        self.cookielog.warn('Creating directories from the User Cookie template.')
        """
        This function uses the username given by our FLASK framework
        and creates a new directory system for the active user using
        our  new_user cookiecutter template.
        """
        if cookie_jar:
            self.cookie_jar = cookie_jar

        # This is used ONLY when the user registers in flask
        # TODO-ROB:  Create the cookiecutter.json file

        # extra_context overrides user and default configs
        cookiecutter(str(self.Recipes.user_cookie), no_input=True, extra_context={
            "user_name": self.user}, output_dir=str(self.cookie_jar))

        # Change user permissions with flask later (this is for testing
        # purposes
        os.chmod(str(self.cookie_jar / Path(self.user)), mode=0o777)
        self.cookielog.info('Directories have been created for the user, %s. ✔' % self.user)

    def bake_the_project(self, cookie_jar=None):
        self.cookielog.warn('Creating directories from the Project Cookie template.')
        """
        :return: A new project inside the user's
        project directory.
        """
        if cookie_jar:
            self.cookie_jar = cookie_jar
        # Add the project
        if self.project:
            no_input = True
            e_c = {"project_name": self.project}
            project_log_message = "(%s)" % self.project
        else:
            no_input = False
            e_c = None
            project_log_message = "that has been named with user input"

        if not self.basic_project:
            self.cookielog.warn('A project linked to a user/repository is being created.')
            cookiecutter(str(self.Recipes.project_cookie), extra_context=e_c, no_input=no_input,
                         output_dir=str(self.cookie_jar))
            # Logging
            if self.user:
                self.cookielog.info('Directories have been created for %s\'s project %s. ✔' % (self.user, project_log_message))
            else:
                self.cookielog.info('Directories have been created for %s.' % project_log_message)
        else:
            self.cookielog.warn('A basic standalone project is being created.')
            cookiecutter(str(self.Recipes.basic_project_cookie), extra_context=e_c, no_input=no_input,
                         output_dir=str(self.cookie_jar))
            self.cookielog.info('Directories have been created for a standalone project %s. ✔' % project_log_message)
        os.chmod(str(self.cookie_jar / Path(self.project)), mode=0o777)

    def bake_the_db_repo(self, db_config_file, db_path, cookie_jar=None, archive_flag=False, delete=False):
        # TODO-ROB:  Work work this in with the database management class.
        """
        :return: A new database inside the users database directory
        """
        if cookie_jar:
            self.cookie_jar = cookie_jar
        # TODO-ROB:  Rework this for new archive function.
        # if archive_flag:
        #     archive_list = archive(database_path=db_path, archive_path=self.cookie_jar, config_file=db_config_file, delete_flag=delete)
        #     for arch in archive_list:
        #         self.cookielog.info("An archive has been created at %s." % arch)
        # else:
        #     if self.db_repo:
        #         no_input = True
        #         e_c = {"db_name": self.db_repo}
        #     else:
        #         no_input = False
        #         e_c = None

            cookiecutter(str(self.Recipes.db_cookie), extra_context=e_c, no_input=no_input, output_dir=str(self.cookie_jar))
            self.cookielog.info("Directories have been created for a database repository %s." %
                                str((self.cookie_jar / Path(self.db_repo))))
            os.chmod(str(self.cookie_jar / Path(self.db_repo)), mode=0o777)
            #
            # for db_key, db_value in db_config_dict["Database_Config"].items():
            #     if db_value:
            #         pass
            # TODO-ROB:  Use db_value system with database management configuration.

    def bake_the_website(self, host, port, website_path, cookie_jar=None):
        self.cookielog.warn('Creating directories from the Website Cookie template.')
        """
        Create a website using the new_website cookie.

        After creating the directory structure, the run_script function
        from cookiecutter finds the hooks folder which contains a
        post-cookiecutter-template-generation bash script.  The bash script
        sets up the proper dependencies and environment variables for the website,
        and runs the website on the specified host and port

        :return: Runs the website.
        """
        if cookie_jar:
            self.cookie_jar = cookie_jar
        # TODO-ROB:  Add heavy logging here
        e_c = {"website_name": self.website,
               "website_path": os.path.join(str(website_path), ''),
               "website_host": host,
               "website_port": port}
        cookiecutter(str(self.Recipes.website_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.cookie_jar))
        os.chmod(str(self.cookie_jar / Path(self.website)), mode=0o777)
        # Get the absolute path to the script that starts the flask server
        script_path = website_path / \
                      Path('hooks') / Path('post_gen_project.sh')
        #scripts_file_path = find_hook('post_gen_project.sh', hooks_dir=str(script_path))
        # TODO-ROB add screening to the bash script for flask run -h -p
        run_script(script_path=str(script_path), cwd=str(website_path))
        self.cookielog.info('Directories have been created for the Flask Web Server, %s. ✔' % self.website)
        self.cookielog.warn('The %s Flask Server should now be running on http://%s:%s' % (self.website, host, port))

    def bake_the_research(self, research_type, research, cookie_jar=None):
        self.cookielog.warn('Creating directories from the Research Cookie template.')
        """
        :param new_app (bool):  Flag for auto generating an app that
         goes with the research target.
        :return:  Adds new directories in the current project labeled
        with the proper names.
        """
        if cookie_jar:
            self.cookie_jar = cookie_jar

        e_c = {"research_type": research_type,
               "research_name": research}
        cookiecutter(str(self.Recipes.research_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.cookie_jar))
        os.chmod(str(self.cookie_jar / Path(research_type)), mode=0o777)
        # script_path = self.project_cookie / Path('hooks') / Path('post_gen_project.py')
        # run_script(script_path, )
        self.cookielog.info('Directories have been created for the %s research project, %s. ✔' % (research_type, research))

    def bake_the_app(self, app, cookie_jar=None):
        self.cookielog.warn('Creating directories from the App Cookie template.')
        """
        Create an app.
        """
        if cookie_jar:
            self.cookie_jar = cookie_jar
        e_c = {"app_name": app}
        cookiecutter(str(self.Recipes.app_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.cookie_jar))
        os.chmod(str(self.cookie_jar), mode=0o777)
        self.cookielog.info("Directories have been created for an R-Shiny app, %s. ✔" % app)
