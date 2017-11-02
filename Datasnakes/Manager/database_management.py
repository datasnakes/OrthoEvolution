import os
from pathlib import Path
from importlib import import_module
from Datasnakes.Manager import ProjectManagement
from Datasnakes.Orthologs.utils import attribute_config
from Datasnakes.Tools.ftp import NcbiFTPClient
from Datasnakes.Tools.logit import LogIt
from Datasnakes.Manager.BioSQL import biosql


class DatabaseManagement(object):

    def __init__(self, project, email, driver, project_path=None, proj_mana=ProjectManagement, biosql_mana=biosql, **kwargs):
        self.dbmanalog = LogIt().default(logname="DatabaseManagement", logfile=None)
        self.config_options = {
            # TODO-ROB:  GI deprecation.  Add some type of window-masker config here.
            # "GI_config": self.get_gi_lists,
            "Blast_config": self.download_blast_database,
            "Taxonomy_config": self.download_taxonomy_database,
            "GenBank_config": {
                "download_taxonomy_database": self.download_taxonomy_database,
                "download_refseq_release_files": self.download_refseq_release_files,
                "upload_refseq_release_files": self.upload_refseq_release_files
                               }
                               }
        self.project = project
        self.email = email
        self.driver = driver
        self.database_dict = {}
        self.ncbiftp = NcbiFTPClient(email=self.email)
        # TODO-ROB:  Configure this differently somehow
        self.biosql = biosql_mana
        self.proj_mana = proj_mana

        # Configuration of class attributes for Project Management.
        if proj_mana:
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
                self.database_dict[db_config_type] = [db_config_method, db_configuration]

# TODO-ROB:  GI deprecation.  Add something for window-masker.
    # def get_gi_lists(self):
    #     print()

    def download_blast_database(self, database_name="refseq_rna", database_path=None):
        # <path>/<user or basic_project>/databases/NCBI/blast/db/<database_name>
        db_path = self.ncbi_db_repo / Path('blast') / Path('db')
        if database_path:
            db_path = Path(database_path)
        self.ncbiftp.getblastdb(database_name=database_name, download_path=self.blast_db)
        # TODO-ROB Add email or slack notifications
        self.dbmanalog.critical("Please set the BLAST environment variables in your .bash_profile!!")
        self.dbmanalog.info("The appropriate environment variable is \'BLASTDB=%s\'." % str(db_path))
        self.dbmanalog.critical("Please set the BLAST environment variables in your .bash_profile!!")
        # TODO-ROB:  set up environment variables.  Also add CLI setup

    def download_taxonomy_database(self, db_type, sub_path):
        """
        This method gets the remote data and updates the local databases for ETE3, BioSQL, and PhyloDB taxonomy
        databases.  Most significant is the "biosql" and "phylodb" types.  The biosql databases use NCBI's taxonomy
        database along with the biosql schema.  And the phylodb databases use ITIS's taxonomy database along with the
        biosql schema.

        :param db_type:  The type of database.  ("ete3", "biosql", or "phylodb")
        :param dest_name:  The name of the new database.
        :param dest_path:  The location where the new database should go.
        :return:  An updated taxonomy database.
        """
        db_type = str(db_type).lower()
        if db_type == 'ete3':
            # DEFAULT_TAXADB = os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')
            ete3 = import_module("ete3")
            ete3.NCBITaxa.update_taxonomy_database()
        elif db_type == 'biosql':
            # Loads data from NCBI via ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy
            if self.driver.lower() == "sqlite3":
                ncbi_db = self.biosql.SQLiteBioSQL(proj_mana=self.proj_mana)
                ncbi_db.copy_template_database(sub_path=sub_path)
                return ncbi_db

            elif self.driver.lower() == "mysql":
                db_path = self.ncbi_db_repo / Path('pub') / Path('taxonomy')
                ncbi_db = self.biosql.MySQLBioSQL()
                return ncbi_db

        elif db_type == 'phylodb':
            # Loads data from ITIS via http://www.itis.gov/downloads/
            print('biosql_repo')

    # TODO-ROB:  Update ncbiftp class syntax to reflect NCBI's ftp site
    def download_refseq_release_files(self, collection_subset, seqtype, filetype):
        db_path = self.ncbi_db_repo / Path('refseq') / Path('release') / Path(collection_subset)
        ncbiftp = self.ncbiftp.getrefseqrelease(database_name=collection_subset, seqtype=seqtype, filetype=filetype, download_path=db_path)
        return ncbiftp

    def upload_refseq_release_files(self, collection_subset, seqtype, filetype, upload_list, extension=".gbk.db"):
        db_path = self.ncbi_refseq_release / Path(collection_subset)
        ncbi_db = self.download_taxonomy_database(db_type="biosql", sub_path="/refseq/release/%s" % collection_subset)
        ncbi_db.upload_path = db_path
        ncbi_db.upload_list = upload_list
        ncbi_db.upload_files(seqtype=seqtype, filetype=filetype)
        pass

    def get_project_genbank_database(self):
        print()

