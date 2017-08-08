from Manager.utils.repo_mana import RepoMana as RM
import os
from pathlib import Path
from cookiecutter.main import cookiecutter

class WebMana(RM):

    def __init__(self, repo, website, host='0.0.0.0', port='5252',
                 home=os.getcwd(), new_website=False, create_admin=False):
        super().__init__(repo=repo, home=home)
        self.website = website
        self.web_host = host
        self.web_port = port
        self.website_path = self.flask / Path(self.website)

        self.website_scripts = self.website_path / Path(self.website)
        self.website_public = self.website_scripts / Path('public')
        self.website_user = self.website_scripts / Path('user')

        if new_website is True:
            self.create_website()

    def create_website(self):
        # TODO-ROB Add heavy logging here
        e_c = {"website_name": self.website,
               "website_path": os.path.join(str(self.website_path), ''),
               "website_host": self.web_host,
               "website_port": self.web_port}
        cookiecutter(str(self.website_cookie), no_input=True,
                     extra_context=e_c, output_dir=self.flask)
        # Get the absolute path to the script that starts the flask server
        script_path = self.website_path / \
            Path('hooks') / Path('post_gen_project.sh')
        #scripts_file_path = find_hook('post_gen_project.sh', hooks_dir=str(script_path))
        run_script(script_path=str(script_path), cwd=str(self.website_path))
