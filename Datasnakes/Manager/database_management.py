#from Datasnakes.Tools.ftp import FTP2DB
import os
from pathlib import Path

from Datasnakes.Manager import ProjectManagement


class DatabaseManagement(object):

    def __init__(self, project, project_path=None, proj_mana=ProjectManagement, **kwargs):

        self.project = project

        if not isinstance(proj_mana, ProjectManagement):
            if project_path:
                self.project_path = Path(project_path) / Path(self.project)
            else:
                self.project_path = Path(os.getcwd()) / Path(self.project)
                Path.mkdir(self.project_path, parents=True, exist_ok=True)
                print('project_path=%s' % self.project_path)
                self.removed_pm_config(kwargs)
        else:
            setattr(proj_mana, 'project', project)
            for key, value in proj_mana.__dict__.items():
                setattr(self, str(key), str(value))
            print('project_path=%s' % self.project_path)

    def removed_pm_config(self, kwargs):
        print()

    def get_gi_lists(self):
        print()

    def get_blast_database(self):
        print()

    def get_taxonomy_database(self):
        print('ete3')
        print('ncbi')
        print('biosql')

    def get_genbank_database(self):
        print()

    def get_project_genbank_database(self):
        print()

