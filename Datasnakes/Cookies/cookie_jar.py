import os
import pkg_resources
from cookiecutter.main import cookiecutter
from pathlib import Path
from Datasnakes import Cookies


class CookieRecipes(object):
    def __init__(self):
        self.CookieJar = Path(pkg_resources.resource_filename(Cookies.__name__, ''))
        self.repo_cookie = self.CookieJar / Path('new_repository')
        self.user_cookie = self.CookieJar / Path('new_user')
        self.project_cookie = self.CookieJar / Path('new_project')
        self.research_cookie = self.CookieJar / Path('new_research')
        self.app_cookie = self.CookieJar / Path('new_app')
        self.db_cookie = self.CookieJar / Path('new_database')
        self.website_cookie = self.CookieJar / Path('new_website')


class Oven(object):

    def __init__(self, repo=None, home=os.getcwd(), cookies=CookieRecipes()):
        self.file_home = home
        # Below are the PyPi path strings
        #    The first group is to access the cookiecutter templates
        self.repo = repo
        self.CookieJar = cookies

    def bake_the_repo(self):
            print('Creating directories from repository cookie.')
            # print(self.__class__.__name__)
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
            cookiecutter(str(self.CookieJar.repo_cookie), no_input=no_input,
                         extra_context=e_c, output_dir=str(self.file_home))
            os.chmod(str(self.file_home / Path(self.repo)), mode=0o777)
            print('Directories have been created. ✔')
