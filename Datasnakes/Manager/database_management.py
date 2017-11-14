from pathlib import Path
import yaml
import requests
from importlib import import_module
from Datasnakes.Manager.management import ProjectManagement
from Datasnakes.Orthologs.utils import attribute_config
from Datasnakes.Cookies.utils import archive
from Datasnakes.Tools.ftp import NcbiFTPClient
from Datasnakes.Tools.logit import LogIt
from Datasnakes.Manager.BioSQL import biosql


class BaseDatabaseManagement(object):

    def __init__(self, project, email, driver, project_path=None, proj_mana=ProjectManagement, **kwargs):
        self.dbmanalog = LogIt().default(logname="DatabaseManagement", logfile=None)
        self.config_options = {
            "WindowMasker_config": self.download_windowmasker_files,
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
        self.biosql = biosql
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

    def download_windowmasker_files(self, taxonomy_ids, database_path=None):
        """Download the WindowMasker files."""
        db_path = self.ncbi_db_repo / Path('blast') / Path('windowmasker_files')
        if database_path:
            db_path = Path(database_path)

        # XXX This will work when we merge.
        self.ncbiftp.getwindowmaskerfiles(taxonomy_ids=taxonomy_ids, download_path=str(db_path))

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

    def download_ncbi_taxonomy_dump_files(self, db_path, url="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"):

        db_abs_path = Path(db_path) / Path('taxdump.tar.gz')
        r = requests.get(url, allow_redirects=True)
        with open(str(db_abs_path), 'wb') as taxdump:
            taxdump.write(r.content)

    # TODO-ROB:  Update ncbiftp class syntax to reflect NCBI's ftp site
    def download_refseq_release_files(self, collection_subset, seqtype, filetype):
        """Download NCBI Refseq Release files."""
        db_path = self.ncbi_db_repo / Path('refseq') / Path('release') / Path(collection_subset)
        ncbiftp = self.ncbiftp.getrefseqrelease(database_name=collection_subset, seqtype=seqtype, filetype=filetype, download_path=db_path)
        return ncbiftp

    def upload_refseq_release_files(self, collection_subset, seqtype, filetype, upload_list, extension=".gbk.db"):
        db_path = self.ncbi_refseq_release / Path(collection_subset)
        # Get a BioSQL database
        ncbi_db = self.download_taxonomy_database(db_type="biosql", sub_path="/refseq/release/%s" % collection_subset)
        ncbi_db.upload_path = db_path
        ncbi_db.upload_list = upload_list
        ncbi_db.upload_files(seqtype=seqtype, filetype=filetype)
        pass

    def get_project_genbank_database(self):
        """"""
        print()


class DatabaseManagement(BaseDatabaseManagement):
    # TODO-ROB: Figure this out for the case of a user.  Because there doesn't necessarily have to be a project
    def __init__(self, config_file, proj_mana=ProjectManagement):
        kw ={}
        db_config_strategy = {}
        with open(config_file) as cf:
            db_config = yaml.load(cf)
            # Get the parameters for the Base class
            for key, value in db_config.items():
                if isinstance(value, dict):
                    db_config_strategy[key] = value
                else:
                    kw[key] = value

        super().__init__(proj_mana=proj_mana, **kw)
        self.strategy_dispatcher = {}
        self.strategy_config = {}
        self.configure_flag = None
        self.archive_flag = None
        self.delete_flag = None
        self.config_file = config_file

    def get_strategy_dispatcher(self, db_config_strategy):
        strategy_dispatcher = {}
        strategy_config = {}
        for strategy, strategy_kwargs in db_config_strategy.items():
            if db_config_strategy[strategy] == "Full":
                strategy_dispatcher, strategy_config = self.full(**strategy_kwargs)
            elif strategy == "NCBI":
                sd, sc = self.ncbi(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "NCBI_blast":
                sd, sc = self.ncbi_blast(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "NCBI_blast_db":
                sd, sc = self.ncbi_blast_db(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "NCBI_blast_windowmaskerfiles":
                sd, sc = self.ncbi_blast_windowmasker_files(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "NCBI_pub_taxonomy":
                sd, sc = self.ncbi_pub_taxonomy(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "NCBI_refseq_release":
                sd, sc = self.ncbi_refseq_release(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "ITIS":
                sd, sc = self.itis(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "ITIS_taxonomy":
                sd, sc = self.itis_taxonomy(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)

        return strategy_dispatcher, strategy_config

    def projects(self, **kwargs):
        print(self)
        return {}, {}

    def full(self, NCBI, ITIS, Projects=None,
             configure_flag=None, archive_flag=None, delete_flag=None):
        if configure_flag:
            NCBI["configure_flag"] = configure_flag
            ITIS["configure_flag"] = configure_flag
        if archive_flag:
            NCBI["archive_flag"] = archive_flag
            ITIS["archive_flag"] = archive_flag
        if delete_flag:
            NCBI["delete_flag"] = delete_flag
            ITIS["delete_flag"] = delete_flag

        full_dispatcher = {}
        full_config = {}
        # Configure NCBI
        ncbi_dispatcher, ncbi_config = self.ncbi(**NCBI)
        # Configure ITIS
        itis_dispatcher, itis_config = self.itis(**ITIS)
        # Configure projects

        # Create Full dispatcher
        full_dispatcher.update(ncbi_dispatcher)
        full_dispatcher.update(itis_dispatcher)
        # Create Full config
        full_config.update(ncbi_config)
        full_config.update(itis_config)
        if Projects:
            projects_dispatcher, projects_config = self.projects(**Projects)
            full_dispatcher.update(projects_dispatcher)
            full_config.update(projects_config)
        # returns dict of config_dicts, dict of dispatcher_functions
        return full_dispatcher, full_config

    def ncbi(self, NCBI_blast, NCBI_pub_taxonomy, NCBI_refseq_release, configure_flag=True, archive_flag=True,
             delete_flag=False, database_path=None, archive_path=None):
        ncbi_dispatcher = {}
        ncbi_config = {}
        # If the flags are set in a top level part of the hierarchy, then everything below follows
        if configure_flag:
            NCBI_blast["configure_flag"] = configure_flag
            NCBI_pub_taxonomy["configure_flag"] = configure_flag
            NCBI_refseq_release["configure_flag"] = configure_flag
        if archive_flag:
            NCBI_blast["delete_flag"] = None
            NCBI_pub_taxonomy["delete_flag"] = None
            NCBI_refseq_release["delete_flag"] = None
            NCBI_blast["archive_flag"] = None
            NCBI_pub_taxonomy["archive_flag"] = None
            NCBI_refseq_release["archive_flag"] = None
            ncbi_dispatcher = {"NCBI": []}
            ncbi_dispatcher["NCBI"].append(archive)
            ncbi_config["NCBI"] = {"NCBI": []}
            ncbi_config["NCBI"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI",
                "delete_flag": delete_flag
            })
        else:
            NCBI_blast["delete_flag"] = delete_flag
            NCBI_pub_taxonomy["delete_flag"] = delete_flag
            NCBI_refseq_release["delete_flag"] = delete_flag

        # Configure blast
        nb_dispatcher, nb_config = self.ncbi_blast(**NCBI_blast)
        # Configure pub/taxonomy
        npt_dispatcher, npt_config = self.ncbi_pub_taxonomy(**NCBI_pub_taxonomy)
        # Configure refseq/release
        nrr_dispatcher, nrr_config = self.ncbi_refseq_release(**NCBI_refseq_release)

        # Create NCBI dispatcher
        ncbi_dispatcher.update(nb_dispatcher)
        ncbi_dispatcher.update(npt_dispatcher)
        ncbi_dispatcher.update(nrr_dispatcher)
        # Create NCBI config
        ncbi_config.update(nb_config)
        ncbi_config.update(npt_config)
        ncbi_config.update(nrr_config)

        return ncbi_dispatcher, ncbi_config

    def ncbi_blast(self, NCBI_blast_db, NCBI_blast_windowmasker_files,
                   configure_flag=None, archive_flag=None, delete_flag=None, database_path=None, archive_path=None):
        ncbi_blast_dispatcher = {}
        ncbi_blast_config = {}
        if configure_flag:
            NCBI_blast_db["configure_flag"] = configure_flag
            NCBI_blast_windowmasker_files["configure_flag"] = configure_flag
        if archive_flag:
            NCBI_blast_db["delete_flag"] = None
            NCBI_blast_windowmasker_files["delete_flag"] = None
            NCBI_blast_db["archive_flag"] = None
            NCBI_blast_windowmasker_files["archive_flag"] = None
            ncbi_blast_dispatcher = {"NCBI_blast": []}
            ncbi_blast_dispatcher["NCBI_blast"].append(archive)
            ncbi_blast_config = {"NCBI_blast": []}
            ncbi_blast_config["NCBI_blast"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_blast",
                "delete_flag": delete_flag
            })
        else:
            NCBI_blast_db["delete_flag"] = delete_flag
            NCBI_blast_windowmasker_files["delete_flag"] = delete_flag

        nbd_dispatcher, nbd_config = self.ncbi_blast_db(**NCBI_blast_db)
        nbw_dispatcher, nbd_config = self.ncbi_blast_windowmasker_files(**NCBI_blast_windowmasker_files)

        ncbi_blast_dispatcher.update(nbd_dispatcher)
        ncbi_blast_dispatcher.update(nbw_dispatcher)

        ncbi_blast_config.update(nbd_config)
        ncbi_blast_config.update(nbw_dispatcher)

        return ncbi_blast_dispatcher, ncbi_blast_config

    def ncbi_blast_db(self, configure_flag=None, archive_flag=None, delete_flag=None, archive_path=None, database_path=None):
        # Set up default parameter values.
        nbd_dispatcher = {"NCBI_blast_db": []}
        nbd_config = {"NCBI_blast_db": []}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        # Archive.  If necessary, then delete.
        if archive_flag:
            nbd_dispatcher["NCBI_blast_db"].append(archive)
            nbd_config["NCBI_blast_db"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_blast",
                "delete_flag": delete_flag
            })
        # Configure
        if configure_flag:
            # Download blast files
            nbd_dispatcher["NCBI_blast_db"].append(self.download_blast_database)
            nbd_config["NCBI_blast_db"].append({
                "database_path": database_path
            })
        return nbd_dispatcher, nbd_config

    def ncbi_blast_windowmasker_files(self, taxonomy_ids, configure_flag=None, archive_flag=None, delete_flag=None,
                                      archive_path=None, database_path=None):
        nbw_dispatcher = {"NCBI_blast_windowmasker_files": []}
        nbw_config = {"NCBI_blast_windowmasker_files": []}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        if archive_flag:
            nbw_dispatcher["NCBI_blast_windowmasker_files"].append(archive)
            nbw_config["NCBI_blast_windowmasker_files"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_blast_windowmasker_files",
                "delete_flag": delete_flag
            })
        if configure_flag:
            nbw_dispatcher["NCBI_blast_windowmasker_files"].append(self.download_windowmasker_files)
            nbw_config["NCBI_blast_windowmasker_files"].append({
                "taxonomy_ids": taxonomy_ids,
                "database_path": database_path
            })
        return nbw_dispatcher, nbw_config

    def ncbi_pub_taxonomy(self, configure_flag=None, archive_flag=None, delete_flag=None, archive_path=None, database_path=None):
        # TODO-ROB:  Add a ftp download of the correct taxdump file for biosql stuff in the self.dl-tax-db method
        npt_dispatcher = {"NCBI_pub_taxonomy": []}
        npt_config = {"NCBI_pub_taxonomy": []}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        if archive_flag:
            npt_dispatcher["NCBI_pub_taxonomy"].append(archive)
            npt_config["NCBI_pub_taxonomy"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_pub_taxonomy",
                "delete_flag": delete_flag
            })
        # Configure
        if configure_flag:
            # Download pub/taxonomy files
            npt_dispatcher["NCBI_blast_db"].append(self.download_ncbi_taxonomy_dump_files)
            npt_config["NCBI_blast_db"].append({
                "db_path": database_path
            })
        return npt_dispatcher, npt_config

    def ncbi_refseq_release(self, configure_flag=None, archive_flag=None, delete_flag=None, upload_flag=None, archive_path=None,
                            database_path=None, collection_subset=None, seqtype=None, filetype=None, upload_list=None):
        nrr_dispatcher = {"NCBI_refseq_release": []}
        nrr_config = {"NCBI_refseq_release": []}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        if archive_flag:
            nrr_dispatcher["NCBI_refseq_release"].append(archive)
            nrr_config["NCBI_refseq_release"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_refseq_release",
                "delete_flag": delete_flag
            })
        if configure_flag:
            nrr_dispatcher["NCBI_refseq_release"].append(self.download_refseq_release_files)
            nrr_config["NCBI_refseq_release"].append({
                "collection_subset": collection_subset,
                "seqtype": seqtype,
                "filetype": filetype
            })
        if upload_flag:
            nrr_dispatcher["NCBI_refseq_release"].append(self.upload_refseq_release_files)
            nrr_config["NCBI_refseq_release"].append({
                "collection_subset": collection_subset,
                "seqtype": seqtype,
                "filetype": filetype,
                "upload_list": upload_list
            })
        return nrr_dispatcher, nrr_config

    def itis(self, ITIS_taxonomy, configure_flag=None, archive_flag=None, delete_flag=None, database_path=None, archive_path=None):
        itis_dispatcher = {}
        itis_config = {}
        if configure_flag:
            ITIS_taxonomy["configure_flag"] = configure_flag
        if archive_flag:
            ITIS_taxonomy["delete_flag"] = None
            ITIS_taxonomy["archive_flag"] = None
            itis_dispatcher = {"ITIS": []}
            itis_dispatcher["ITIS"].append(archive)
            itis_config = {"ITIS": []}
            itis_config["ITIS"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "ITIS",
                "delete_flag": delete_flag
            })
        else:
            ITIS_taxonomy["delete_flag"] = delete_flag
        it_dispatcher, it_config = self.itis_taxonomy(**ITIS_taxonomy)
        itis_dispatcher.update(it_dispatcher)
        itis_config.update(it_config)

        return itis_dispatcher, itis_config

    def itis_taxonomy(self, configure_flag=None, archive_flag=None, delete_flag=None, **kwargs):
        return {}, {}
