import os
from pathlib import Path
from importlib import import_module
from Datasnakes.Manager import ProjectManagement
from Datasnakes.Orthologs.utils import attribute_config
from Datasnakes.Tools.ftp import NcbiFTPClient
from Datasnakes.Tools.logit import LogIt


class DatabaseManagement(object):

    def __init__(self, project, email, project_path=None, proj_mana=ProjectManagement, **kwargs):
        self.dbmanalog = LogIt().default(logname="DatabaseManagement", logfile=None)
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

    def get_blast_database(self, database_name="refseq_rna", database_path=None):
        # <path>/<user or basic_project>/databases/NCBI/blast/db/<database_name>
        db_path = self.ncbi_db_repo / Path('blast') / Path('db')
        if database_path:
            db_path = Path(database_path)
        self.ncbiftp.getblastdb(database_name=database_name, download_path=db_path)
        # TODO-ROB Add email or slack notifications
        self.dbmanalog.critical("Please set the BLAST environment variables in your .bash_profile!!")
        self.dbmanalog.info("The appropriate environment variable is \'BLASTDB=%s\'." % str(db_path))
        self.dbmanalog.critical("Please set the BLAST environment variables in your .bash_profile!!")
        # TODO-ROB:  set up environment variables.  Also add CLI setup

    def get_taxonomy_database(self, db_type, db_name=None, dest_name=None, dest_path=None, driver=None):
        db_type = str(db_type).lower()
        if db_type == 'ete3':
            # DEFAULT_TAXADB = os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')
            ete3 = import_module("ete3")
            ete3.NCBITaxa.update_taxonomy_database()
        elif db_type == 'biosql':
            # Loads data from NCBI via ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy
            biosql = import_module("Datasnakes.Manager.BioSQL.biosql")
            if driver.lower() == "sqlite":
                db_path = self.ncbi_db_repo / Path('pub') / Path('taxonomy')
                ncbi_db = biosql.SQLiteBioSQL(database_name=db_name, database_path=db_path)
                ncbi_db.copy_template_database(db_path)

            elif driver.lower() == "mysql":
                db_path = self.ncbi_db_repo / Path('pub') / Path('taxonomy')
                ncbi_db = biosql.MySQLBioSQL()

        elif db_type == 'phylodb':
            # Loads data from ITIS via http://www.itis.gov/downloads/
            print('biosql_repo')

    def get_genbank_database(self):
        print()

    def get_project_genbank_database(self):
        print()

