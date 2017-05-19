import os
from pathlib import Path

from cookiecutter.main import cookiecutter

from archive.utils.user_mana import UserMana as UM


class ProjMana(UM):

    def __init__(self, repo, user, project, research=None, research_type=None, app=None, home=os.getcwd(),
                 new_project=False, new_research=False, new_app=False):
        """
        :param repo: 
        :param user: 
        :param project: 
        :param research: 
        :param app: 
        :param home: 
        :param new_project: 
        :param project_type: 
        :param new_research: 
        :param new_app: 
        """
        super().__init__(repo=repo, user=user, project=project, home=home,
                         new_project=new_project)

        # TODO-ROB The init for this is wrong.  THis was done fro the new_research cookie.
        # TODO-ROB Go back to the drawing board for the public/private/other choices.  (FLASK forms)
        self.project = project
        self.research = research
        self.research_type = research_type

        # Project/Research Directories
        self.research_path = self.project_path / Path(research_type) / Path(research)
        self.data = self.research_path / Path('data')
        self.raw_data = self.research_path / Path('raw_data')
        self.project_web = self.research_path / Path('web')
        # TODO-ROB:  THis is just a draft.  Rework to use public/private/other
        if app:
            self.app = app
            self.app_path = self.project_web / Path(app)
        if new_research is True:
            self.create_research(new_app)

    def create_research(self, new_app=False):
        e_c = {"research_type": self.research_type,
               "research_name": self.research}
        cookiecutter(self.research_cookie, no_input=True, extra_context=e_c, output_dir=self.research_path)
        if new_app is True:
            self.create_app()

    def create_app(self):
        e_c = {"app_name": self.app}
        cookiecutter(self.app_cookie, no_input=True, extra_context=e_c, output_dir=self.app_path)
