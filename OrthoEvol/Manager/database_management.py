# Standard Library
import os
import random
import string
import urllib.request
import tarfile
import yaml
from collections import OrderedDict
from importlib import import_module
from pathlib import Path
from pkg_resources import resource_filename
import subprocess as sp
# OrthoEvol
from OrthoEvol import OrthoEvolDeprecationWarning
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.Tools.ftp import NcbiFTPClient
from OrthoEvol.utilities import FullUtilities
from OrthoEvol.Manager.biosql import biosql
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Orthologs.Blast.comparative_genetics import BaseComparativeGenetics
from OrthoEvol.Manager.config import templates


class BaseDatabaseManagement(object):

    def __init__(self, email, driver, project=None, project_path=None, proj_mana=None, blast=False, ftp_flag=True):
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
        :param ftp_flag: A flag used if FTP connection is available or not.
        :type ftp_flag:  bool.
        """

        # Initialize Utilities
        self.db_mana_utils = FullUtilities()
        self.db_mana_log = LogIt().default(logname="Database-Management", logfile=None)
        self.project = project
        self.email = email
        self.driver = driver
        self.database_dict = {}
        self.ftp_flag = ftp_flag

        if ftp_flag:
            try:
                self.ncbiftp = NcbiFTPClient(email=self.email)
            except:
                self.ftp_flag = False
                self.db_mana_log.warning("This system doesn't allow FTP usage.")

        self.biosql = biosql
        self.proj_mana = proj_mana

        # Configuration of class attributes for Project Management.
        if proj_mana:
            add_self = self.db_mana_utils.attribute_config(self, composer=proj_mana, checker=ProjectManagement, project=project, project_path=project_path)
            for var, attr in add_self.__dict__.items():
                setattr(self, var, attr)
            self.database_path = self.user_db
            if blast:
                self.gene_data = BaseComparativeGenetics(project=self.project, project_path=project_path, proj_mana=proj_mana, copy_from_package=True, MAF='MAFV3.3_short.csv')
        else:
            self.database_path = Path(project_path) / Path(project) / Path("databases")

    def download_windowmasker_files(self, taxonomy_ids):
        """Download the WindowMasker files used in the BLAST database.

        :param taxonomy_ids:  Taxonomy ids for the organisms of interest.
        :type taxonomy_ids:  list.
        """
        OrthoEvolDeprecationWarning("Windowmasker files are no longer used by the most current blastx command line "
                                    "utilities.  You can now use taxon ids directly.")
        # <path>/<user or basic_project>/databases/NCBI/blast/db/<database_name>
        dl_path = Path(self.database_path) / Path("NCBI") / Path('blast') / Path('windowmasker_files')
        self.ncbiftp.getwindowmaskerfiles(taxonomy_ids=taxonomy_ids, download_path=str(dl_path))

    def download_blast_database(self, database_name="refseq_rna", v5=True, set_blastdb=True):
        """Download the blast database files for using NCBI's BLAST+ command line.

        For other types of blast data, please see the NCBIREADME.md file.

        :param database_name:  A string that represents a pattern in the files of interest, defaults to "refseq_rna"
        :type database_name:  str, optional
        :param v5: A flag that determines which version of blastdb to use, defaults to True
        :type v5: bool, optional
        :param set_blastdb:  A flag that determines whether the BLASTDB environment
                            variable is automatically set.
        :type set_blastdb:  bool, optional
        """
        # <path>/<user or basic_project>/databases/NCBI/blast/db/<database_name>
        dl_path = Path(self.database_path) / Path("NCBI") / Path("blast") / Path("db")
        
        if v5:
            dl_path = Path(dl_path) / Path("v5")

        # Download the preformatted blast database.
        self.ncbiftp.getblastdb(database_name=database_name, v5=v5,
                                download_path=str(dl_path))

        env_vars = dict(os.environ).keys()
        if set_blastdb or ("BLASTDB" not in env_vars):
            # See if .bash_profile or .profile exists
            bash_prof = Path("~/.bash_profile").expanduser().absolute()
            sh_prof = Path("~/.profile").expanduser().absolute()
            if not bash_prof.exists():
                if not sh_prof.exists():
                    bash_prof.touch(mode=0o700)
                    set_prof = bash_prof
                else:
                    set_prof = sh_prof
            else:
                set_prof = bash_prof
            self.db_mana_log.warning("Setting the PATH in %s" % str(set_prof))
            # Use the set .*profile to append to PATH
            with open(str(set_prof), 'r') as prof:
                _ = prof.read()
                bas_prof_export = "export PATH=\"%s:$PATH\"" % str(dl_path)
                if bas_prof_export not in _:
                    with open(str(set_prof), "a+") as b_prof:
                        b_prof.write("export PATH=\"%s:$PATH\"" % str(dl_path))
                    cmd = ["source %s" % str(set_prof)]
                    stdout = self.db_mana_utils.system_cmd(cmd=cmd,
                                                           stdout=sp.PIPE,
                                                           stderr=sp.STDOUT,
                                                           shell=True)
        else:
            self.db_mana_log.critical("Please set the BLAST environment variables in your .bash_profile!!")
            self.db_mana_log.info("The appropriate environment variable is \'BLASTDB=%s\'." % str(dl_path))
            self.db_mana_log.critical("Please set the BLAST environment variables in your .bash_profile!!")

    def download_ete3_taxonomy_database(self):
        """Update ETE3's taxonomy database with ETE3's API."""
        # DEFAULT_TAXADB = os.path.join(os.environ.get('HOME', '/'), '.etetoolkit', 'taxa.sqlite')
        ete3 = import_module("ete3")
        ncbi_taxon_dump_path = self.database_path / Path("NCBI") / Path('pub') / Path('taxonomy')
        ncbi = ete3.NCBITaxa(dbfile=str(ncbi_taxon_dump_path / 'ete3_taxa.sqlite'))
        ncbi.update_taxonomy_database()

    def copy_biosql_taxonomy_template(self, destination, database_name):
        """
        Copy a BioSQL template database loaded with NCBI's taxonomy data.
        :param destination:  Where the template will be copied to.
        :type destination:  str.
        :param database_name:  The name of the copied database.
        :type database_name:  str.
        """
        if self.driver.lower() == "sqlite3":
            ncbi_db = self.biosql.SQLiteBioSQL(database_name=database_name,
                                               proj_mana=self.proj_mana)
            ncbi_db.copy_template_database(destination=destination)
            return ncbi_db

        elif self.driver.lower() == "mysql":
            # db_path = self.ncbi_db_repo / Path('pub') / Path('taxonomy')
            # ncbi_db = self.biosql.MySQLBioSQL()
            # return ncbi_db
            pass

    def create_biosql_taxonomy_template(self):
        """Creates a template database by uploading SQLite schema and NCBI taxonomy."""
        if self.driver.lower() == "sqlite3":
            ncbi_db = self.biosql.SQLiteBioSQL(proj_mana=self.proj_mana)
            ncbi_db.create_template_database()
        return ncbi_db

    def download_itis_taxonomy_tables(self):
        # Loads data from ITIS via http://www.itis.gov/downloads/
        # Use this along with BioSQL's phyloDB
        pass

    def download_ncbi_taxonomy_dump_files(self, url='''ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'''):
        """Download and extract the NCBI taxonomy dump files via a GET request.

        :param url: A ftp link to the NCBI taxdump*.tar.gz file of interest.
        :type url: str.
        """
        # TODO: Add this to NCBI FTP client
        dl_path = Path(self.database_path) / Path("NCBI") / Path('pub') / Path('taxonomy')
        dl_abs_path = dl_path / Path('taxdump.tar.gz')
        r = urllib.request.urlopen(url)
        data = r.read()
        with open(str(dl_abs_path), 'wb') as taxdump:
            taxdump.write(data)
        with tarfile.open(str(dl_abs_path)) as tar:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(tar, dl_path)
        os.remove(dl_abs_path)

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
        :return: A list of files to download from NCBI via FTP.
        :rtype:  list.
        """
        db_path = self.database_path / Path('NCBI') / Path('refseq') / Path('release') / Path(collection_subset)
        db_path.mkdir(parents=True, exist_ok=True)
        # TODO: If database exists and is same size, use the existing database.
        self.ncbiftp.getrefseqrelease(collection_subset=collection_subset, seqtype=seqtype,
                                      seqformat=seqformat, download_path=db_path)
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
        :param upload_list:  A list of files to upload.
        :type upload_list:  list.
        :param database_name:  The name of the database to create.  The default name is usually best.
        :type database_name:  str.
        :param add_to_default:  A string to add to the default name.
        :type add_to_default:  str.
        """

        if database_name:
            db_name = database_name
        else:
            if add_to_default:
                add_to_default = "_%s" % add_to_default
            else:
                add_to_default = ""
            db_name = "{}_{}{}.{}.db".format(collection_subset, seqtype, add_to_default, seqformat)
        db_path = self.database_path / Path("NCBI") / Path("refseq") / Path("release") / Path(collection_subset)
        # Get a BioSQL database
        ncbi_db = self.copy_biosql_taxonomy_template(destination=db_path, database_name=db_name)
        ncbi_db.upload_files(seqtype=seqtype, filetype='genbank', upload_path=db_path, upload_list=upload_list)

    def get_project_genbank_database(self):
        """"""
        pass


class DatabaseManagement(BaseDatabaseManagement):

    def __init__(self, config_file, proj_mana=ProjectManagement, **kwargs):
        """
        This class creates higher level functionality for configuring various databases.  It uses a YAML configuration
        file to parse various database setup "strategies" and dispatches or serializes the setup functions in the
        proper order.  The strategies range from very general ("Full", "NCBI", etc.) to very specific ("NCBI_blast_db",
        "NCBI_refseq_release", etc.).  Each strategy has 3 separate flags (configure, archive, and delete).  Specific
        strategies may also have other flags or parameters that control their behavior.
        :param config_file:  The path to a YAML configuration file.  See databases.yml in the package config
        folder.
        :type config_file:  str.
        :param proj_mana: A configuration variable for connecting projects.
        :type proj_mana: ProjectManagement.
        :param kwargs:  Key-word arguments.
        :type kwargs:  dict.
        """
        self.db_mana_utils = FullUtilities()
        self.db_config_strategy, kw = self.db_mana_utils.parse_db_config_file(config_file)
        super().__init__(proj_mana=proj_mana, **kw)
        self.strategy_dispatcher = OrderedDict()
        self.strategy_config = OrderedDict()
        self.configure_flag = None
        self.archive_flag = None
        self.delete_flag = None
        self.config_file = config_file

    def get_strategy_dispatcher(self, db_config_strategy):
        """
        Loop through a dictionary of strategies with nested configurations, and return a list of functions, and a
        list of matching key-word arguments (dictionaries).  The functions can then be dispatched using the kwargs.
        Higher level more generalized workflows will have lists nested within dictionaries.

        :param db_config_strategy:  A dictionary usually generated from a YAML configuration file.
        :type db_config_strategy:  dict.
        :return:  A tuple containing 2 objects:  a list of functions, and a list of dictionaries containing kwargs for
        each function.
        :rtype:  tuple.
        """
        strategy_dispatcher = OrderedDict()
        strategy_config = OrderedDict()
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

    def full(self, NCBI, ITIS, Projects=None, configure_flag=None, archive_flag=None, delete_flag=None, project_flag=None, _path=None):
        """
        The most generalized strategy available.  This configures everything.  The 3 primary flags (configure, archive,
        and delete) will be passed down to the more specific strategies, which will inherit these values unless
        expressly overridden.

        :param NCBI:  Keyword arguments for the generalized NCBI strategy.
        :type NCBI:  dict.
        :param ITIS:  Keyword arguments for the generalized ITIS strategy.
        :type ITIS:  dict.
        :param configure_flag:  A generalized flag that is passed to all of the strategies in order to implement their
        configuration process.
        :type configure_flag:  bool.
        :param archive_flag:  A generalized flag that is passed to all of the strategies in order to implement their
        archiving process.
        :type archive_flag:  bool.
        :param delete_flag:  A generalized flag that is passed to all of the strategies in order to implement their
        deletion process.
        :type delete_flag:  bool.
        :param _path:
        :type _path:
        :return:  A tuple containing 2 objects:  a list of function (NCBI and ITIS), and a list of dictionaries
        containing kwargs for each function.  In the future Projects will also be  returned.
        :rtype:  tuple.
        """
        if configure_flag:
            NCBI["configure_flag"] = configure_flag
            ITIS["configure_flag"] = configure_flag
        if archive_flag:
            NCBI["archive_flag"] = archive_flag
            ITIS["archive_flag"] = archive_flag
        if delete_flag:
            NCBI["delete_flag"] = delete_flag
            ITIS["delete_flag"] = delete_flag

        full_dispatcher = OrderedDict()
        full_config = OrderedDict()
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

    def NCBI(self, NCBI_blast, NCBI_pub_taxonomy, NCBI_refseq_release, configure_flag=None, archive_flag=None,
             delete_flag=None, database_path=None, archive_path=None, _path=None):
        """
        A strategy that implements all of the databases relevant to NCBI.

        :param NCBI_blast:  Keyword arguments for the generalized NCBI_blast strategy.
        :type NCBI_blast:  dict.
        :param NCBI_pub_taxonomy:  Keyword arguments for the specific NCBI_pub_taxonomy strategy.
        :type NCBI_pub_taxonomy:  dict.
        :param NCBI_refseq_release:  Keyword arguments for the specific NCBI_refseq_release strategy.
        :type NCBI_refseq_release:  dict.
        :param configure_flag:  A generalized flag that is passed to the NCBI strategies in order to implement their
        configuration process.
        :type configure_flag:  bool.
        :param archive_flag:  A generalized flag that is passed to to the NCBI strategies in order to implement their
        archiving process.
        :type archive_flag:  bool.
        :param delete_flag:  A generalized flag that is passed to the NCBI strategies in order to implement their
        deletion process.
        :type delete_flag:  bool.
        :param database_path:  User supplied relative path to the databases.
        :type database_path:   str.
        :param archive_path:  User supplied relative path to the archived databases.
        :type archive_path:   str.
        :return:  A tuple containing 2 objects:  a list of function (NCBI_blast, NCBI_pub_taxonomy, and
        NCBI_refseq_release), and a list of dictionaries containing kwargs for each function.
        :rtype:  tuple.
        """
        ncbi_dispatcher = OrderedDict({"NCBI": []})
        ncbi_config = OrderedDict({"NCBI": []})
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
        """
        A strategy that implements all of the data relevant to NCBI's blast databases.

        :param NCBI_blast_db:  Keyword arguments for the generalized NCBI_blast_db strategy.
        :type NCBI_blast_db:   dict.
        :param NCBI_blast_windowmasker_files:  Keyword arguments for the generalized NCBI_blast_windowmasker_files
         strategy.
        :type NCBI_blast_windowmasker_files:  dict.
        :param configure_flag:  A flag that is passed to the NCBI strategies in order to implement their
        configuration process.
        :type configure_flag:  bool.
        :param archive_flag:  A flag that is passed to to the NCBI strategies in order to implement their
        archiving process.
        :type archive_flag:  bool.
        :param delete_flag:  A flag that is passed to the NCBI strategies in order to implement their
        deletion process.
        :type delete_flag:  bool.
        :param database_path:  User supplied relative path to the databases.
        :type database_path:   str.
        :param archive_path:  User supplied relative path to the archived databases.
        :type archive_path:   str.
        :return:  A tuple containing 2 objects:  a list of function (NCBI_blast_db, and NCBI_blast_windowmasker_files),
        and a list of dictionaries containing kwargs for each function.
        :rtype:  tuple.
        """
        ncbi_blast_dispatcher = OrderedDict()
        ncbi_blast_config = OrderedDict()
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
        """
        A strategy that implements NCBI's blast database that's used with the blast+ command line utilities.

        :param configure_flag:  A flag for configuring the blast database for NCBI's blast+ tool.
        :type configure_flag:  bool.
        :param archive_flag:  A flag for archiving the blast database for NCBI's blast+ tool.
        :type archive_flag:  bool.
        :param delete_flag:  A flag for deleting the blast database for NCBI's blast+ tool.
        :type delete_flag:  bool.
        :param database_path:  User supplied relative path to the databases.
        :type database_path:   str.
        :param archive_path:  User supplied relative path to the archived databases.
        :type archive_path:   str.
        :return:  A tuple containing 2 objects:  a list of functions for dealing with NCBI's blast database,
        and a list of dictionaries containing kwargs for each function.
        :rtype:  tuple.
        """
        # Set up default parameter values.
        nbd_dispatcher = OrderedDict({"NCBI_blast_db": []})
        nbd_config = OrderedDict({"NCBI_blast_db": []})
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
        """
        A strategy that sets up windowmasker files used with the blast+ command line utilities.

        :param configure_flag:  A flag for configuring the windowmasker files for blast+.
        :type configure_flag:  bool.
        :param archive_flag:  A flag for archiving the windowmasker files for blast+.
        :type archive_flag:  bool.
        :param delete_flag:  A flag for deleting the windowmasker files for blast+.
        :type delete_flag:  bool.
        :param database_path:  User supplied relative path to the databases.
        :type database_path:   str.
        :param archive_path:  User supplied relative path to the archived databases.
        :type archive_path:   str.
        :return:  A tuple containing 2 objects:  a list of functions for dealing with windowmasker files,
        and a list of dictionaries containing kwargs for each function.
        :rtype:  tuple.
        """
        nbw_dispatcher = OrderedDict({"NCBI_blast_windowmasker_files": []})
        nbw_config = OrderedDict({"NCBI_blast_windowmasker_files": []})
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
        """
        A strategy that sets up windowmasker files used with the blast+ command line utilities.

        :param configure_flag:  A flag that updates NCBI's taxonomy dump files.
        :type configure_flag:  bool.
        :param archive_flag:  A flag that archives NCBI's taxonomy dump files.
        :type archive_flag:  bool.
        :param delete_flag:  A flag that deletes NCBI's taxonomy dump files.
        :type delete_flag:  bool.
        :param database_path:  User supplied relative path to the databases.
        :type database_path:   str.
        :param archive_path:  User supplied relative path to the archived databases.
        :type archive_path:   str.
        :return:  A tuple containing 2 objects:  a list of functions for downloading NCBI's taxdum.tar.gz,
        and a list of dictionaries containing kwargs for each function.
        :rtype:  tuple.
        """
        npt_dispatcher = OrderedDict({"NCBI_pub_taxonomy": []})
        npt_config = OrderedDict({"NCBI_pub_taxonomy": []})
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
                            upload_number=8, _path=None, activate=None, template_flag=None, download_flag=None,
                            pbs_dict=None):
        """
        This is the most complicated specific strategy.  It downloads the refseq release files of choice (gbff),
        extracts the data, splits a list of files into {upload_number} lists, and uploads those file lists to
        {upload_number} BioSQL databases.  The resulting databases can be used to access sequence data with accession
        numbers.  Due to the long upload time (>24 hours) this class currently only uses PBS to break the uploading up
        into {upload_number} processes.

        :param configure_flag:  A flag that downloads refseq release files from NCBI.
        :type configure_flag:  bool.
        :param archive_flag:  A flag that archives refseq release files from NCBI.
        :type archive_flag:  bool.
        :param delete_flag:  A flag that deletes refseq release files from NCBI.
        :type delete_flag:  bool.
        :param upload_flag:  A flag that uploads refseq release files to BioSQL database(s).
        :type upload_flag:  bool.
        :param database_path:  User supplied relative path to the databases.
        :type database_path:   str.
        :param archive_path:  User supplied relative path to the archived databases.
        :type archive_path:   str.
        :param collection_subset: The collection subset of interest.
        :type collection_subset: str.
        :param seqtype: The type of sequence (rna, protein, genomic).
        :type seqtype: str.
        :param seqformat: The format of the sequence file (usually 'gbff' for GenBank Flat File).
        :type seqformat: str.
        :param file_list:  A list of files to upload.
        :type file_list:  list.
        :param upload_number:  The number of databases to upload (defautls to 8)
        :type upload_number:  int.
        :param activate:  Absolute path to a virtual environments activate script.  This is used for PBS scripts.
        :type activate:  str.
        :param template_flag:  A flag that loads a BioSQL database with NCBI taxonomy data.  This takes a very long time.
        :type template_flag:  bool.
        :param download_flag:  A flag that downloads the proper refseq release files from NCBI's ftp site.
        :type download_flag:  bool.
        :return:  A tuple containing 2 objects:  a list of functions for managing a BioSQL database with NCBI refseq
        data, and a list of dictionaries containing kwargs for each function.
        :rtype:  tuple.
        """
        nrr_dispatcher = OrderedDict({"NCBI_refseq_release": OrderedDict({"archive": [], "configure": [], "upload": []})})
        nrr_config = OrderedDict({"NCBI_refseq_release": OrderedDict({"archive": [], "configure": [], "upload": []})})
        _biosql = self.biosql.SQLiteBioSQL(proj_mana=self.proj_mana)
        dl_path = Path(self.database_path) / Path("NCBI") / Path('pub') / Path('taxonomy')
        dmp_file = dl_path / "nodes.dmp"
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
            if download_flag and self.ftp_flag:
                nrr_dispatcher["NCBI_refseq_release"]["configure"].append(self.download_refseq_release_files)
                nrr_config["NCBI_refseq_release"]["configure"].append({
                    "collection_subset": collection_subset,
                    "seqtype": seqtype,
                    "seqformat": seqformat
                })
            if template_flag and (self.ftp_flag or dmp_file.exists()):
                nrr_dispatcher["NCBI_refseq_release"]["configure"].append(self.create_biosql_taxonomy_template)
                nrr_config["NCBI_refseq_release"]["configure"].append({})

        if upload_flag:
            if not _biosql.template_abs_path.is_file() and (self.ftp_flag or dmp_file.exists()):
                nrr_dispatcher["NCBI_refseq_release"]["configure"].append(self.create_biosql_taxonomy_template)
                nrr_config["NCBI_refseq_release"]["configure"].append({})
            elif _biosql.template_abs_path.is_file():
                self.db_mana_log.info("The BioSQL template exists.")
            elif not self.ftp_flag and dmp_file.exists():
                self.db_mana_log.info("The system does not allow FTP for downloading the taxonomic dump files, but"
                                      "they already exist.")
            else:
                self.db_mana_log.error("The BioSQL template doesn't exist and the system does not allow FTP for "
                                       "downloading the taxonomic dump files.  The dump files are also not already"
                                       "available.")

            if upload_number < 8:
                raise ValueError("The upload_number must be greater than 8.  The NCBI refseq release files are too bing"
                                 "for anything less than 8 seperate BioSQL databases.")

            # Get template script variables
            py_shebang = Path(activate.parent).expanduser()
            db_path = self.database_path / Path('NCBI') / Path('refseq') / Path('release') / Path(collection_subset)

            # Read the upload script
            upload_script = resource_filename(templates.__name__, 'upload_rr_pbs.py')
            with open(upload_script, 'r') as upload_script:
                temp_script = upload_script.read()
            rand_str = random.sample(string.ascii_letters + string.digits, 5)
            script_dir = Path(self.user_log, ('upload_rr' + ''.join(rand_str)))
            script_dir.mkdir()
            script_string = temp_script % (py_shebang, file_list, pbs_dict, db_path, upload_number, self.email, str(script_dir / 'upload_config.yml'))

            # Create the master upload script
            with open(str(script_dir / 'master_upload_rr_pbs.py'), 'w') as mus:
                mus.write(script_string)
            os.chmod(str(script_dir / 'master_upload_rr_pbs.py'), mode=0o755)

            # Load the configuration file
            with open(str(self.config_file), 'r') as cfg:
                conf_data = yaml.load(cfg, Loader=yaml.FullLoader)

            # Write to the upload config file
            with open(str(script_dir / 'upload_config.yml'), 'w') as upload_cfg:
                conf_data['Database_config']['ftp_flag'] = False
                conf_data['Database_config']['Full']['NCBI']['NCBI_refseq_release']['upload_flag'] = False
                yaml.dump(conf_data, upload_cfg, default_flow_style=False)
            os.chmod(str(script_dir / 'upload_config.yml'), mode=0o755)

            def _run_upload_script():
                self.db_mana_utils.system_cmd(cmd='%s/master_upload_rr_pbs.py' % str(script_dir), cwd=str(script_dir),
                                              stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)

            nrr_dispatcher["NCBI_refseq_release"]['upload'].append(_run_upload_script)
            nrr_config["NCBI_refseq_release"]['upload'].append({})

        return nrr_dispatcher, nrr_config

    def itis(self, ITIS_taxonomy, configure_flag=None, archive_flag=None, delete_flag=None, database_path=None, archive_path=None, _path=None):
        itis_dispatcher = OrderedDict()
        itis_config = OrderedDict()
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
