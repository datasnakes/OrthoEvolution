import os
from pathlib import Path
from cookiecutter.main import cookiecutter
from Manager.utils.repo_mana import RepoMana as RM


class UserMana(RM):

    def __init__(self, user, home=os.getcwd(), project=None, new_project=False):
        super().__init__(home=home, user=user)
        self.user = user
        self.user_path = self.user_path

        self.user_index = self.user_path / Path('index')
        self.user_log = self.user_path / Path('log')
        self.manuscripts = self.user_path / Path('manuscripts')
        self.other = self.user_path / Path('other')
        self.projects = self.user_path / Path('projects')

        if project:
            self.project = project
            self.project_path = self.user_path / Path(project)
        if new_project is True:
            self.create_project()

    def create_project(self):
        """
        :return: A new project inside the user's
        project directory.
        """
        cookiecutter(self.project_cookie, no_input=True, extra_context={"new_project": self.project},
                     output_dir=self.project_path)

