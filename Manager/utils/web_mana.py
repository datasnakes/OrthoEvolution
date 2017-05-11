from Manager.utils.repo_mana import RepoMana as RM
import os
from pathlib import Path

class WebMana(RM):

    def __init__(self, repo, website, home=os.getcwd(), new_website=False):
        super().__init__(repo=repo, home=home)
        self.website = website
        self.website_path = self.flask / Path(self.website)

        self.website_scripts = self.website_path / Path(self.website)
        self.website_public = self.website_scripts / Path('public')
        self.website_user = self.website_scripts / Path('user')

        if new_website is True:
            self.create_website()

    def create_website(self):
        print('web')