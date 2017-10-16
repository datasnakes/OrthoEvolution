import os
from pathlib import Path

from Datasnakes.Manager import ProjectManagement
from Datasnakes.Orthologs.utils import attribute_config


class DatabaseManagement(object):

    def __init__(self, project, project_path=None, proj_mana=ProjectManagement, **kwargs):
        self.config_options = {
            "GI_config": self.get_gi_lists,
            "Blast_config": self.get_blast_database,
            "Taxonomy_config": self.get_taxonomy_database,
            "GenBank_config": self.get_genbank_database
                               }
        self.project = project

        # Configuration of class attributes.
        add_self = attribute_config(self, composer=proj_mana, checker=ProjectManagement, project=project, project_path=project_path)
        for var, attr in add_self.__dict__.items():
            setattr(self, var, attr)

        # Determine which database to update
        # And then run that script with the configuration.
        for config in self.config_options.keys():
            if config in kwargs.keys():
                db_config_method = self.config_options[config]
                db_config_method(kwargs[config])

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

