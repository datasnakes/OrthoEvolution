import os
from pathlib import Path
from cookiecutter.main import cookiecutter

from Manager.utils.mana import Mana as DM


class RepoMana(DM):

    def __init__(self, repo, home=os.getcwd(), user=None, new_user=False):
        super().__init__(home=home, repo=repo)

        self.docs = self.repo_path / Path('docs')
        self.misc = self.repo_path / Path('misc')
        self.users = self.repo_path / Path('users')

        self.lib = self.repo_path / Path('lib')
        self.archives = self.lib / Path('archives')
        self.databases = self.lib / Path('databases')
        self.repo_index = self.lib / Path('index')

        self.repo_web = self.repo_path / Path('web')
        self.flask = self.repo_web / Path('flask')
        self.repo_shiny = self.repo_web / Path('shiny')
        self.ftp = self.repo_web / Path('ftp')
        self.wasabi = self.repo_web / Path('wasabi')

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
        print('hold')
        pass

    def git_ignore(self):
        """Get the ignored file patterns from the .gitignore file in the repo."""
        with open(self.repo_path + '.gitignore', 'r', newline='') as ignore:
            ignored = ignore.read().splitlines()
        return ignored



# TODO-ROB:  Edit the setup.py file for cookiecutter.
