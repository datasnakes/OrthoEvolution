from cookiecutter.main import cookiecutter
import os
from pathlib import Path
import sqlite3, sqlalchemy, flask_sqlalchemy  # one of these
# TODO-ROB in flask add project registration along with user registration

cookiecutter(Path(os.environ['ACTIVE_USER'])/Path('projects'), no_input=True, extra_context=os.environ['EXTRA_CONTEXT'],
             output_dir=os.environ['ACTIVE_USER']/'projects')
cookiecutter()