import os
from pathlib import Path

from Datasnakes.Manager import ProjectManagement
from Datasnakes.Orthologs.utils import attribute_config
from Datasnakes.Tools.ftp import NcbiFTPClient


class DatabaseManagement(object):

    def __init__(self, project, email, project_path=None, proj_mana=ProjectManagement, **kwargs):
        self.config_options = {
            "GI_config": self.get_gi_lists,
            "Blast_config": self.get_blast_database,
            "Taxonomy_config": self.get_taxonomy_database,
            "GenBank_config": self.get_genbank_database
                               }
        self.project = project
        self.email = email
        self.database_dict = {}
        self.ncbiftp = NcbiFTPClient(email=self.email)

        # Configuration of class attributes.
        add_self = attribute_config(self, composer=proj_mana, checker=ProjectManagement, project=project, project_path=project_path)
        for var, attr in add_self.__dict__.items():
            setattr(self, var, attr)

        # Determine which database to update
        # And then run that script with the configuration.
        for config in self.config_options.keys():
            if config in kwargs.keys():
                db_config_type = config
                db_config_method = self.config_options[config]
                db_configuration = kwargs[config]
                # db_config_method(kwargs[config])
                self.database_dict[db_config_type] = [db_config_method, db_configuration]

    def get_gi_lists(self):
        print()

    def get_blast_database(self, database_name):
        # <path>/<user or basic_project>/databases/NCBI/blast/db/<database_name>
        db_path = self.ncbi_db_repo / Path('blast') / Path('db') / Path(database_name)
        self.ncbiftp.getblastdb(database_name=database_name, download_path=db_path)

    def get_taxonomy_database(self, db_type):
        db_type = str(db_type).lower()
        if db_type == 'ete3':
            print('ete3')
        elif db_type == 'ncbi':
            print('ncbi')
        elif db_type == 'biosql':
            print('biosql')

    def get_genbank_database(self):
        print()

    def get_project_genbank_database(self):
        print()

