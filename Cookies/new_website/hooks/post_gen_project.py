import os
from pathlib import Path
from Manager.utils.repo_mana import RepoMana
#TODO-ROB Create a Website_mana class
Man = RepoMana()
flask_path = Man.flask

os.system('source post_gen_project.sh')
