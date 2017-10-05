import os
import pkg_resources
from cookiecutter.main import cookiecutter
from cookiecutter.prompt import prompt_for_config
from cookiecutter.generate import generate_context
from cookiecutter.hooks import run_script
from pathlib import Path
from Datasnakes import Cookies
from Datasnakes.Tools import LogIt


class CookieRecipes(object):
    def __init__(self):
        self.CookieJar = Path(pkg_resources.resource_filename(Cookies.__name__, ''))
        self.repo_cookie = self.CookieJar / Path('new_repository')
        self.user_cookie = self.CookieJar / Path('new_user')
        self.project_cookie = self.CookieJar / Path('new_project')
        self.basic_project_cookie = self.CookieJar / Path('new_basic_project')
        self.research_cookie = self.CookieJar / Path('new_research')
        self.app_cookie = self.CookieJar / Path('new_app')
        self.db_cookie = self.CookieJar / Path('new_database')
        self.website_cookie = self.CookieJar / Path('new_website')


class Oven(object):

    def __init__(self, repo=None, user=None, project=None, basic_project=False, databases=None, website=None, output_dir=os.getcwd(), cookies=CookieRecipes()):
        self.cookielog = LogIt().default(logname="Cookies", logfile=None)
        self.cookie_jar = output_dir
        # Below are the PyPi path strings
        #    The first group is to access the cookiecutter templates
        self.repo = repo
        self.user = user
        self.project = project
        self.basic_project = basic_project
        self.databases = databases
        self.website = website
        self.Ingredients = cookies

    def bake_the_repo(self, cookie_jar=None):
            self.cookielog.warn('Creating directories from the Repository Cookie template.')
            # print(self.__class__.__name__)
            """This function creates a new repository.  If a repository name
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
            cookiecutter(str(self.Ingredients.repo_cookie), no_input=no_input,
                         extra_context=e_c, output_dir=str(self.cookie_jar))
            os.chmod(str(self.cookie_jar / Path(self.repo)), mode=0o777)
            self.cookielog.info('Repository directories have been created. ✔')

    def bake_the_user(self, cookie_jar=None):
        self.cookielog.warn('Creating directories from the User Cookie template.')
        """This function uses the username given by our FLASK framework
        and creates a new directory system for the active user using
        our  new_user cookiecutter template.
        """
        if cookie_jar:
            self.cookie_jar = cookie_jar

        # This is used ONLY when the user registers in flask
        # TODO-ROB:  Create the cookiecutter.json file

        # extra_context overrides user and default configs
        cookiecutter(str(self.Ingredients.user_cookie), no_input=True, extra_context={
            "user_name": self.user}, output_dir=str(self.cookie_jar))

        # Change user permissions with flask later (this is for testing
        # purposes
        os.chmod(str(self.cookie_jar / Path(self.user)), mode=0o777)
        self.cookielog.info('Directories have been created for the user, %s. ✔' % self.user)

    def bake_the_project(self, cookie_jar=None):
        self.cookielog.warn('Creating directories from the Project Cookie template.')
        # print(self.__class__.__name__)
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
            cookiecutter(str(self.Ingredients.project_cookie), extra_context=e_c, no_input=no_input,
                         output_dir=str(self.cookie_jar))
            # Logging
            if self.user:
                self.cookielog.info('Directories have been created for %s\'s project %s. ✔' % (self.user, project_log_message))
            else:
                self.cookielog.info('Directories have been created for %s.' % project_log_message)
        else:
            self.cookielog.warn('A basic standalone project is being created.')
            cookiecutter(str(self.Ingredients.basic_project_cookie), extra_context=e_c, no_input=no_input,
                         output_dir=str(self.cookie_jar))
            self.cookielog.info('Directories have been created for a standalone project %s. ✔' % project_log_message)
        os.chmod(str(self.cookie_jar / Path(self.project)), mode=0o777)

    def bake_the_db_repo(self, user_db, db_path_dict=None, ncbi_db_repo=None):
        self.cookielog.warn('Creating directories from the Database Cookie.')
        """
        :return: A new database inside the users database directory
        """
        # TODO-ROB:  FIx this.  output_dir needs to take self.cookie_jar
        # TODO-ROB:  There is a better way to accomplish this
        if db_path_dict:
            for db, path in db_path_dict.items():
                e_c = {"db_name": db}
                cookiecutter(str(self.Ingredients.db_cookie), extra_context=e_c, no_input=True, output_dir=str(user_db))
                os.chmod(str(user_db / Path(db)), mode=0o777)
                self.cookielog.info('Directories have been created for the database, %s. ✔' % db)

        else:
            db_num = int(input("How many NCBI databases do you need to create?"))
            for db in range(1, db_num + 1):
                # Manually set up cookiecutter prompting
                context_file = str(self.Ingredients.db_cookie / Path('cookiecutter.json'))
                e_c = prompt_for_config(context=generate_context(context_file=context_file))
                # Create the cookiecutter repo with no input, and add extra content from manual prompts
                cookiecutter(str(self.Ingredients.db_cookie), output_dir=str(ncbi_db_repo), extra_context=e_c, no_input=True)
                # Use cookiecutter_dict from manual prompts to change the user permissions.
                os.chmod(str(ncbi_db_repo / Path(e_c['db_name'])), mode=0o777)
                self.cookielog.info('Directories have been created for the database, %s. ✔' % e_c['db_name'])

    def bake_the_website(self, host, port, website_path, cookie_jar=None):
        self.cookielog.warn('Creating directories from the Website Cookie template.')
        """Create a website using the new_website cookie.

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
        cookiecutter(str(self.Ingredients.website_cookie), no_input=True,
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
        cookiecutter(str(self.Ingredients.research_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.cookie_jar))
        os.chmod(str(self.cookie_jar / Path(research_type)), mode=0o777)
        # script_path = self.project_cookie / Path('hooks') / Path('post_gen_project.py')
        # run_script(script_path, )
        self.cookielog.info('Directories have been created for the %s research project, %s. ✔' % (research_type, research))

    def bake_the_app(self, app, cookie_jar=None):
        self.cookielog.warn('Creating directories from the App Cookie template.')
        """Create an app."""
        if cookie_jar:
            self.cookie_jar = cookie_jar
        e_c = {"app_name": app}
        cookiecutter(str(self.Ingredients.app_cookie), no_input=True,
                     extra_context=e_c, output_dir=str(self.cookie_jar))
        os.chmod(str(self.cookie_jar), mode=0o777)
        self.cookielog.info("Directories have been created for an R-Shiny app, %s. ✔" % app)
