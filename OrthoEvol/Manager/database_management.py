# Standard Library
import os
import urllib.request
import tarfile
from importlib import import_module
from pathlib import Path
# OrthoEvol
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Tools.ftp import NcbiFTPClient
from OrthoEvol.utilities import FullUtilities
from OrthoEvol.Manager.BioSQL import biosql
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Orthologs.Blast.comparative_genetics import BaseComparativeGenetics
from OrthoEvol.Tools.sge import SGEJob


class BaseDatabaseManagement(object):

    def __init__(self, email, driver, project=None, project_path=None, proj_mana=ProjectManagement, ftp_flag=True):
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

        # Initialize Utilities
        self.db_mana_utils = FullUtilities()
        self.db_mana_log = LogIt().default(logname="DatabaseManagement", logfile=None)
        self.project = project
        self.email = email
        self.driver = driver
        self.database_dict = {}
        if ftp_flag == True:
            self.ncbiftp = NcbiFTPClient(email=self.email)
        self.biosql = biosql
        self.proj_mana = proj_mana

        # Configuration of class attributes for Project Management.
        if proj_mana:
            add_self = self.db_mana_utils.attribute_config(self, composer=proj_mana, checker=ProjectManagement, project=project, project_path=project_path)
            for var, attr in add_self.__dict__.items():
                setattr(self, var, attr)
            self.database_path = self.user_db
            self.gene_data = BaseComparativeGenetics(project=self.project, project_path=project_path, proj_mana=proj_mana, copy_from_package=True, MAF='MAFV3.3_short.csv')
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

    def create_biosql_taxonomy_database(self, destination, database_name):
        """Create a BioSQL database template loaded with NCBI's taxonomy data."""
        if self.driver.lower() == "sqlite3":
            ncbi_db = self.biosql.SQLiteBioSQL(database_name=database_name, proj_mana=self.proj_mana)
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

    def download_ncbi_taxonomy_dump_files(self, url='''ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'''):
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
        r = urllib.request.urlopen(url)
        data = r.read()
        print('taxdump.tar.gz')
        with open(str(dl_abs_path), 'wb') as taxdump:
            taxdump.write(data)
        with tarfile.open(str(dl_abs_path)) as tar:
            tar.extractall(dl_path)
        os.remove(dl_abs_path)

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
        db_path.mkdir(parents=True, exist_ok=True)
        ncbiftp = self.ncbiftp.getrefseqrelease(collection_subset=collection_subset, seqtype=seqtype, seqformat=seqformat, download_path=db_path)
        return self.ncbiftp.files2download

    def upload_refseq_release_files(self, collection_subset, seqtype, seqformat, upload_list=None, database_name=None, add_to_default=None):
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

        if database_name:
            db_name = database_name
        else:
            if add_to_default:
                add_to_default = "_%s" % add_to_default
            db_name = "{}_{}{}.{}.db".format(collection_subset, seqtype, add_to_default, seqformat)
        db_path = self.database_path / Path("NCBI") / Path("refseq") / Path("release") / Path(collection_subset)
        # Get a BioSQL database
        ncbi_db = self.create_biosql_taxonomy_database(destination=db_path, database_name=db_name)
        ncbi_db.upload_files(seqtype=seqtype, filetype='genbank', upload_path=db_path, upload_list=upload_list)

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
    def __init__(self, config_file, proj_mana=ProjectManagement, **kwargs):
        self.db_mana_utils = FullUtilities()
        self.db_config_strategy, kw = self.db_mana_utils.parse_db_config_file(config_file)
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
            if strategy == "Full":
                strategy_dispatcher, strategy_config = self.full(**strategy_kwargs)
            elif strategy == "NCBI":
                sd, sc = self.NCBI(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "NCBI_blast":
                sd, sc = self.NCBI_blast(**strategy_kwargs)
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
                sd, sc = self.NCBI_pub_taxonomy(**strategy_kwargs)
                strategy_dispatcher.update(sd)
                strategy_config.update(sc)
            elif strategy == "NCBI_refseq_release":
                sd, sc = self.NCBI_refseq_release(**strategy_kwargs)
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

    def PROJECTS(self, **kwargs):
        print(self)
        return {}, {}

    def full(self, NCBI, ITIS, Projects=None,
             configure_flag=None, archive_flag=None, delete_flag=None, project_flag=None, _path=None):
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
        ncbi_dispatcher, ncbi_config = self.NCBI(**NCBI)
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
            projects_dispatcher, projects_config = self.PROJECTS(**Projects)
            full_dispatcher.update(projects_dispatcher)
            full_config.update(projects_config)
        # returns dict of config_dicts, dict of dispatcher_functions
        return full_dispatcher, full_config

    def NCBI(self, NCBI_blast, NCBI_pub_taxonomy, NCBI_refseq_release, configure_flag=True, archive_flag=True,
             delete_flag=False, database_path=None, archive_path=None, _path=None):
        ncbi_dispatcher = {"NCBI": []}
        ncbi_config = {"NCBI": []}
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
            ncbi_dispatcher["NCBI"].append(self.db_mana_utils.archive)
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
        nb_dispatcher, nb_config = self.NCBI_blast(**NCBI_blast)
        # Configure pub/taxonomy
        npt_dispatcher, npt_config = self.NCBI_pub_taxonomy(**NCBI_pub_taxonomy)
        # Configure refseq/release
        nrr_dispatcher, nrr_config = self.NCBI_refseq_release(**NCBI_refseq_release)

        # Create NCBI dispatcher
        ncbi_dispatcher.update(nb_dispatcher)
        ncbi_dispatcher.update(npt_dispatcher)
        ncbi_dispatcher.update(nrr_dispatcher)
        # Create NCBI config
        ncbi_config.update(nb_config)
        ncbi_config.update(npt_config)
        ncbi_config.update(nrr_config)

        return ncbi_dispatcher, ncbi_config

    def NCBI_blast(self, NCBI_blast_db, NCBI_blast_windowmasker_files,
                   configure_flag=None, archive_flag=None, delete_flag=None, database_path=None, archive_path=None, _path=None):
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
            ncbi_blast_dispatcher["NCBI_blast"].append(self.db_mana_utils.archive)
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
        nbw_dispatcher, nbw_config = self.ncbi_blast_windowmasker_files(**NCBI_blast_windowmasker_files)

        ncbi_blast_dispatcher.update(nbd_dispatcher)
        ncbi_blast_dispatcher.update(nbw_dispatcher)

        ncbi_blast_config.update(nbd_config)
        ncbi_blast_config.update(nbw_config)

        return ncbi_blast_dispatcher, ncbi_blast_config

    def ncbi_blast_db(self, configure_flag=None, archive_flag=None, delete_flag=None, archive_path=None, database_path=None, _path=None):
        # Set up default parameter values.
        nbd_dispatcher = {"NCBI_blast_db": []}
        nbd_config = {"NCBI_blast_db": []}
        if not archive_path:
            archive_path = str(self.user_archive)
        # Archive.  If necessary, then delete.
        if archive_flag:
            nbd_dispatcher["NCBI_blast_db"].append(self.db_mana_utils.archive)
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
                "database_name": 'refseq_rna'
            })
        return nbd_dispatcher, nbd_config

    def ncbi_blast_windowmasker_files(self, taxonomy_ids, configure_flag=None, archive_flag=None, delete_flag=None,
                                      archive_path=None, database_path=None, _path=None):
        nbw_dispatcher = {"NCBI_blast_windowmasker_files": []}
        nbw_config = {"NCBI_blast_windowmasker_files": []}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        if archive_flag:
            nbw_dispatcher["NCBI_blast_windowmasker_files"].append(self.db_mana_utils.archive)
            nbw_config["NCBI_blast_windowmasker_files"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_blast_windowmasker_files",
                "delete_flag": delete_flag
            })
        if configure_flag:
            nbw_dispatcher["NCBI_blast_windowmasker_files"].append(self.download_windowmasker_files)
            nbw_config["NCBI_blast_windowmasker_files"].append({
                "taxonomy_ids": self.gene_data.taxon_ids
            })
        return nbw_dispatcher, nbw_config

    def NCBI_pub_taxonomy(self, configure_flag=None, archive_flag=None, delete_flag=None, archive_path=None, database_path=None, _path=None):
        # TODO-ROB:  Add a ftp download of the correct taxdump file for biosql stuff in the self.dl-tax-db method
        npt_dispatcher = {"NCBI_pub_taxonomy": []}
        npt_config = {"NCBI_pub_taxonomy": []}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        if archive_flag:
            npt_dispatcher["NCBI_pub_taxonomy"].append(self.db_mana_utils.archive)
            npt_config["NCBI_pub_taxonomy"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_pub_taxonomy",
                "delete_flag": delete_flag
            })
        # Configure
        if configure_flag:
            # Download pub/taxonomy files
            npt_dispatcher["NCBI_pub_taxonomy"].append(self.download_ncbi_taxonomy_dump_files)
            npt_config["NCBI_pub_taxonomy"].append({})
        return npt_dispatcher, npt_config

    def NCBI_refseq_release(self, configure_flag=None, archive_flag=None, delete_flag=None, upload_flag=None, archive_path=None,
                            database_path=None, collection_subset=None, seqtype=None, seqformat=None, file_list=None,
                            upload_number=8, _path=None, activate=None):
        nrr_dispatcher = {"NCBI_refseq_release": {"archive": [], "configure": [], "upload": []}}
        nrr_config = {"NCBI_refseq_release": {"archive": [], "configure": [], "upload": []}}
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        if archive_flag:
            nrr_dispatcher["NCBI_refseq_release"]["archive"].append(self.db_mana_utils.archive)
            nrr_config["NCBI_refseq_release"]["archive"].append({
                "database_path": database_path,
                "archive_path": archive_path,
                "option": "NCBI_refseq_release",
                "delete_flag": delete_flag
            })
        if configure_flag:
            nrr_dispatcher["NCBI_refseq_release"]["configure"].append(self.download_refseq_release_files)
            nrr_config["NCBI_refseq_release"]["configure"].append({
                "collection_subset": collection_subset,
                "seqtype": seqtype,
                "seqformat": seqformat
            })
        if upload_flag:
            if upload_number < 8:
                raise ValueError("The upload_number must be greater than 8.  The NCBI refseq release files are too bing"
                                 "for anything less than 8 seperate BioSQL databases.")

            # Create a list of lists with an index corresponding to the upload number
            if file_list is None:
                db_path = self.database_path / Path('NCBI') / Path('refseq') / Path('release') / Path(collection_subset)
                file_list = os.listdir(str(db_path))
                file_list = [x for x in file_list if x.endswith(str(seqformat))]
            sub_upload_size = len(file_list) // upload_number
            sub_upload_lists = [file_list[x:x + 100] for x in range(0, len(file_list), sub_upload_size)]
            if (len(file_list) % upload_number) != 0:
                upload_number = upload_number + 1
            add_to_default = 0
            for sub_list in sub_upload_lists:
                add_to_default += 1
                nrr_dispatcher["NCBI_refseq_release"]["upload"].append(self.db_mana_utils.refseq_jobber)
                code_dict_string = str({
                    "collection_subset": collection_subset,
                    "seqtype": seqtype,
                    "seqformat": seqformat,
                    "upload_list": sub_list,
                    "add_to_default": add_to_default
                })
                # Create a Python script for this in the package
                sge_code_string = \
                "from OrthoEvol.Manager.management import ProjectManagement\n" \
                "from OrthoEvol.Manager.database_dispatcher import DatabaseDispatcher\n" \
                "from OrthoEvol.Manager.config import yml\n" \
                "from pkg_resources import resource_filename\n" \
                "import yaml\n" \
                "pm_config_file = resource_filename(yml.__name__, \"config_template_existing.yml\")\n" \
                "with open(pm_config_file, \'r\') as f:\n" \
                "   pm_config = yaml.safe_load(f)\n" \
                "pm = ProjectManagement(**pm_config[\"Management_config\"])\n" \
                "code_dict_string = %s\n" \
                "R_R = DatabaseDispatcher(config_file=\"%s\", proj_mana=pm, upload_refseq_release=True, **code_dict_string)\n" % \
                    (code_dict_string, self.config_file)
                nrr_config["NCBI_refseq_release"]["upload"].append({
                    "code": sge_code_string,
                    "base_jobname": "upload_rr_%s",
                    "email_address": self.email,
                    "id": add_to_default,
                    "activate": activate})

        return nrr_dispatcher, nrr_config

    def itis(self, ITIS_taxonomy, configure_flag=None, archive_flag=None, delete_flag=None, database_path=None, archive_path=None, _path=None):
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
            itis_dispatcher["ITIS"].append(self.db_mana_utils.archive)
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


class DatabaseDispatcher(DatabaseManagement):

    def __init__(self, config_file, **kwargs):
        super().__init__(config_file=config_file)
        self.dispatcher, self.configuration = self.get_strategy_dispatcher(db_config_strategy=kwargs)
        if len(self.dispatcher.keys()) == 1:
            self.strategy = list(self.dispatcher.keys())[0]
        else:
            raise ValueError

        self.archive_disp = self.dispatcher[self.strategy]["archive"]
        self.configure_disp = self.dispatcher[self.strategy]["configure"]
        self.upload_disp = self.dispatcher[self.strategy]["upload"]

        self.archive_config = self.configuration[self.strategy]["archive"]
        self.configure_config = self.configuration[self.strategy]["configure"]
        self.upload_config = self.configuration[self.strategy]["configure"]

    def dispatch_the_uploader(self):
        self.upload_disp(**self.upload_config)