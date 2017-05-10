import os
from pathlib import Path
from cookiecutter.main import cookiecutter
from Manager.utils.user_mana import UserMana as UM


class ProjMana(UM):

    def __init__(self, user, project, home=os.getcwd(), research=None, app=None, new_research=False, new_app=False):
        """

        :param user:
        :param project:
        :param home:
        :param research:
        :param new_research:
        """
        super().__init__(user=user, project=project, home=home)
        # TODO-ROB The init for this is wrong.  THis was done fro the new_research cookie.
        # TODO-ROB Go back to the drawing board for the public/private/other choices.  (FLASK forms)
        self.user = user
        self.user_path = self.user_path
        self.project = project
        self.project_path = self.project_path

        self.data = self.project_path / Path('data')
        self.raw_data = self.project_path / Path('raw_data')
        self.project_web = self.project_path / Path('web')
        # TODO-ROB:  THis is just a draft.  Rework to use public/private/other
        if research:
            self.research = research
            self.research_path = self.project_path / Path(research)
            if app:
                self.app = app
                self.app_path = self.project_web / Path(app)
        if new_research is True:
            self.create_research()
            if new_app is True:
                self.create_app()

    def create_research(self):

        cookiecutter(self.research_cookie, no_input=True, extra_context={"new_research": self.research},
                     output_dir=self.research_path)

    def create_app(self):

        cookiecutter(self.app_cookie, no_input=True, extra_context={"new_app": self.app}, output_dir=self.app_path)
