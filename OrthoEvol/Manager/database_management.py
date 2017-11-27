"""Database management tools for the OrthoEvol package."""
from pathlib import Path
import yaml
import requests
import tarfile
from importlib import import_module
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Orthologs.utils import attribute_config
from OrthoEvol.Cookies.utils import archive
from OrthoEvol.Tools.ftp import NcbiFTPClient
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Manager.BioSQL import biosql


class BaseDatabaseManagement(object):

    def __init__(self, project, email, driver, project_path=None, proj_mana=ProjectManagement):
        """
        This is the base class for managing various databases.  It provides functionality for downloading and creating
        various databases for your pipeline.  There are functions available for downloading files from NCBI (BLAST,
        windowmasker, taxonomy, refseq release), downloading ITIS taxonomy tables, creating BioSQL databases, and
        uploading refseq release files to BioSQL databases.  This class currently REQUIRES an instance of
        ProjectManagement to be used with the proj_mana parameter.

        :param project: The name of the project.
        :type project: str.
        :param email: The email of the user for using during the FTP.
        :type email: str.
        :param driver: The driver used for creating the BioSQL databases.
        :type driver:  str.
        :param project_path: A path used for standalone/basic project configuration.
        :type project_path: str.
        :param proj_mana: A configuration variable for connecting projects.
        :type proj_mana: ProjectManagement.
        """

        self.db_mana_log = LogIt().default(logname="DatabaseManagement", logfile=None)
        self.project = project
        self.email = email
        self.driver = driver
        self.database_dict = {}
        self.ncbiftp = NcbiFTPClient(email=self.email)
        self.biosql = biosql
        self.proj_mana = proj_mana

        # Configuration of class attributes for Project Management.
        if proj_mana:
            add_self = attribute_config(self, composer=proj_mana, checker=ProjectManagement, project=project, project_path=project_path)
            for var, attr in add_self.__dict__.items():
                setattr(self, var, attr)
            self.database_path = self.user_db
        else:
            self.database_path = Path(project_path)

    def download_windowmasker_files(self, taxonomy_ids):
        """
        Download the WindowMasker files used in the BLAST database.

        :param taxonomy_ids:  Taxonomy ids for the organisms of interest.
        :type taxonomy_ids:  list.
        :return:
        :rtype:
        """
        # <path>/<user or basic_project>/databases/NCBI/blast/db/<database_name>
        dl_path = Path(self.database_path) / Path("NCBI") / Path('blast') / Path('windowmasker_files')
        self.ncbiftp.getwindowmaskerfiles(taxonomy_ids=taxonomy_ids, download_path=str(dl_path))

    def download_blast_database(self, database_name="refseq_rna"):
        """
        Download the blast database files for using NCBI's BLAST+ command line.  For other types of blast data, please
        see the NCBIREADME.md file.

        :param database_name:  A string that represents a pattern in the files of interest.
        :type database_name:  string.
        :return:
        :rtype:
        """
        # <path>/<user or basic_project>/databases/NCBI/blast/db/<database_name>
        dl_path = Path(self.database_path) / Path("NCBI") / Path("blast") / Path("db")

        self.ncbiftp.getblastdb(database_name=database_name, download_path=str(dl_path))

        # TODO-ROB Add email or slack notifications
        self.dbmanalog.critical("Please set the BLAST environment variables in your .bash_profile!!")
        self.dbmanalog.info("The appropriate environment variable is \'BLASTDB=%s\'." % str(dl_path))
        self.dbmanalog.critical("Please set the BLAST environment variables in your .bash_profile!!")
        # TODO-ROB:  set up environment variables.  Also add CLI setup

    def download_ete3_taxonomy_database(self):
        """Update ETE3's taxonomy database with ETE3's API."""
        # DEFAULT_TAXADB = os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')
        ete3 = import_module("ete3")
        ete3.NCBITaxa.update_taxonomy_database()

    def create_biosql_taxonomy_database(self, destination):
        """Create a BioSQL database template loaded with NCBI's taxonomy data."""
        if self.driver.lower() == "sqlite3":
            ncbi_db = self.biosql.SQLiteBioSQL(proj_mana=self.proj_mana)
            ncbi_db.copy_template_database(destination=destination)
            return ncbi_db

        elif self.driver.lower() == "mysql":
            # db_path = self.ncbi_db_repo / Path('pub') / Path('taxonomy')
            # ncbi_db = self.biosql.MySQLBioSQL()
            # return ncbi_db
            pass

    def download_itis_taxonomy_tables(self):
        # Loads data from ITIS via http://www.itis.gov/downloads/
        # Use this along with BioSQL's phyloDB
        pass

    def download_ncbi_taxonomy_dump_files(self, url="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"):
        """
        Download and extract the NCBI taxonomy dump files via a GET request.

        :param url: A ftp link to the NCBI taxdump.* file of interest.
        :type url: str.
        :return:
        :rtype:
        """
        # TODO-ROB:  Maybe add taxdata folder here to reflect BioSQL's load_ncbi_taxonomy.pl script
        dl_path = Path(self.database_path) / Path("NCBI") / Path('pub') / Path('taxonomy')
        dl_abs_path = dl_path / Path('taxdump.tar.gz')
        r = requests.get(url, allow_redirects=True)
        with open(str(dl_abs_path), 'wb') as taxdump:
            taxdump.write(r.content)
        with tarfile.open(str(dl_abs_path)) as tar:
            tar.extractall(dl_path)

    # TODO-ROB:  Update ncbiftp class syntax to reflect NCBI's ftp site
    def download_refseq_release_files(self, collection_subset, seqtype, seqformat):
        """
        Download NCBI Refseq Release files from NCBI.  The collection subtype is a species group
        (e.g. vertebrate_mammalian) of interest.

        :param collection_subset: The collection subset of interest.
        :type collection_subset: str.
        :param seqtype: The type of sequence (RNA, protein, genomic).
        :type seqtype: str.
        :param seqformat: The format of the sequence file (usually 'gbff' for GenBank Flat File).
        :type seqformat: str.
        :return:
        :rtype:
        """
        db_path = self.database_path / Path('NCBI') / Path('refseq') / Path('release') / Path(collection_subset)
        ncbiftp = self.ncbiftp.getrefseqrelease(taxon_group=collection_subset, seqtype=seqtype, seqformat=seqformat, download_path=db_path)
        return ncbiftp

    def upload_refseq_release_files(self, collection_subset, seqtype, seqformat, upload_number=8, upload_list=None, extension=".gbk.db"):
        """
        Upload NCBI's Refseq Release files to a BioSQL database.

        :param collection_subset: The collection subset of interest.
        :type collection_subset: str.
        :param seqtype: The type of sequence (RNA, protein, genomic).
        :type seqtype: str.
        :param seqformat: The format of the sequence file (usually 'gbff' for GenBank Flat File).
        :type seqformat: str.
        :param upload_list:
        :type upload_list:
        :param extension:
        :type extension:
        :return:
        :rtype:
        """
        if upload_number < 8:
            raise ValueError("The upload_number must be greater than 8.  The NCBI refseq release files are too bing"
                             "for anything less than 8 seperate BioSQL databases.")
        from OrthoEvol.Tools.sge import SGEJob

        # TODO-ROB: multiprocessing sge whatever
        # Create a list of lists with an index corresponding to the upload number
        sub_upload_size = len(upload_list) / upload_number
        sub_upload_lists = [upload_list[x:x+100] for x in range(0, len(upload_list), sub_upload_size)]
        upload_job = SGEJob(email_address=self.email)
        for sub_list in sub_upload_lists:
            db_path = self.database_path / Path("NCBI") / Path("refseq") / Path("release") / Path(collection_subset)
            # Get a BioSQL database
            ncbi_db = self.create_biosql_taxonomy_database(destination=db_path)
            ncbi_db.upload_files(seqtype=seqtype, filetype=seqformat, upload_path=db_path, upload_list=upload_list)

    def get_project_genbank_database(self):
        """"""
        print()

    # def download_taxonomy_database(self, db_type, destination):
    #     """
    #     This method gets the remote data and updates the local databases for ETE3, BioSQL, and PhyloDB taxonomy
    #     databases.  Most significant is the "biosql" and "phylodb" types.  The biosql databases use NCBI's taxonomy
    #     database along with the biosql schema.  And the phylodb databases use ITIS's taxonomy database along with the
    #     biosql schema.
    #
    #     :param db_type:  The type of database.  ("ete3", "biosql", or "phylodb")
    #     :param dest_name:  The name of the new database.
    #     :param dest_path:  The location where the new database should go.
    #     :return:  An updated taxonomy database.
    #     """
        # db_type = str(db_type).lower()
        # if db_type == 'ete3':
        #     # DEFAULT_TAXADB = os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')
        #     ete3 = import_module("ete3")
        #     ete3.NCBITaxa.update_taxonomy_database()
        # elif db_type == 'biosql':
        #     # Loads data from NCBI via ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy
        #     if self.driver.lower() == "sqlite3":
        #         ncbi_db = self.biosql.SQLiteBioSQL(proj_mana=self.proj_mana)
        #         ncbi_db.copy_template_database(destination=destination)
        #         return ncbi_db
        #
        #     elif self.driver.lower() == "mysql":
        #         db_path = self.ncbi_db_repo / Path('pub') / Path('taxonomy')
        #         ncbi_db = self.biosql.MySQLBioSQL()
        #         return ncbi_db

        # elif db_type == 'phylodb':
        #     # Loads data from ITIS via http://www.itis.gov/downloads/
        #     print('biosql_repo')


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
                    continue
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
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
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
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
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
                "database_path": database_path
            })
        return npt_dispatcher, npt_config

    def ncbi_refseq_release(self, configure_flag=None, archive_flag=None, delete_flag=None, upload_flag=None, archive_path=None,
                            database_path=None, collection_subset=None, seqtype=None, seqformat=None, upload_list=None,
                            upload_number=8):
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
                "seqformat": seqformat
            })
        if upload_flag:
            if upload_number < 8:
                raise ValueError("The upload_number must be greater than 8.  The NCBI refseq release files are too bing"
                                 "for anything less than 8 seperate BioSQL databases.")
            from OrthoEvol.Tools.sge import SGEJob

            # TODO-ROB: multiprocessing sge whatever
            # Create a list of lists with an index corresponding to the upload number
            sub_upload_size = len(upload_list) / upload_number
            if (len(upload_list) % upload_number) != 0:
                upload_number = upload_number + 1
            sub_upload_lists = [upload_list[x:x+100] for x in range(0, len(upload_list), sub_upload_size)]
            upload_job = SGEJob(email_address=self.email)
            for sub_list in sub_upload_lists:
                self.upload_refseq_release_files(collection_subset=collection_subset, seqtype=seqtype, seqformat=seqformat, upload_list=sub_list)
                nrr_dispatcher["NCBI_refseq_release"].append(self.upload_refseq_release_files)
                nrr_config["NCBI_refseq_release"].append({
                    "collection_subset": collection_subset,
                    "seqtype": seqtype,
                    "seqformat": seqformat,
                    "upload_list": sub_list
                })
        return nrr_dispatcher, nrr_config

    def itis(self, ITIS_taxonomy, configure_flag=None, archive_flag=None, delete_flag=None, database_path=None, archive_path=None):
        itis_dispatcher = {}
        itis_config = {}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
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
