import os
from pathlib import Path

from cookiecutter.main import cookiecutter

from archive.utils.repo_mana import RepoMana as RM


class UserMana(RM):
    #         # TODO-ROB CREATE THESE IN A VIRTUAL ENVIRONMENT FOR EACH USER
    #         # TODO-ROB The virtual environment can be the name of the user
    #         # TODO-ROB When the user logs in, they will activate the virtual environment
    #     # TODO-ROB USE SQL here to see if the user db contains the username
    def __init__(self, repo, user, project=None, home=os.getcwd(), new_user=False, new_project=False):
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
