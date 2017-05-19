import os
from pathlib import Path
from cookiecutter.main import cookiecutter

from Manager.utils import mana
DM = mana.Mana


class RepoMana(DM):

    def __init__(self, repo, user=None, home=os.getcwd(), new_user=False, new_repo=False):
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
        # TODO-ROB:  This is used ONLY when the user registers in flask
        # TODO-ROB:  Create the cookiecutter.json file
        # extra_context overrides user and default configs
        cookiecutter(self.user_cookie, no_input=True, extra_context={"user_name": self.user}, output_dir=self.users)



# TODO-ROB:  Edit the setup.py file for cookiecutter.
