# Standard Library
import os
import subprocess as sp
import tarfile
import urllib.request
from collections import OrderedDict
from collections.abc import Callable, Sequence
from importlib import import_module
from pathlib import Path
from typing import Any

import yaml

# OrthoEvol
from OrthoEvol import OrthoEvolDeprecationWarning
from OrthoEvol.Manager.biosql import biosql
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Orthologs.Blast.comparative_genetics import \
    BaseComparativeGenetics
from OrthoEvol.Tools.ftp import NcbiFTPClient
from OrthoEvol.Tools.logit import LogIt
from OrthoEvol.utilities import FullUtilities


class BaseDatabaseManagement(object):

    def __init__(self, email, driver, project=None, project_path=None, proj_mana=None, blast=False, ftp_flag=True):
        """
        This is the base class for managing various databases.  It provides functionality for downloading and creating
        various databases for your pipeline.  There are functions available for downloading files from NCBI (BLAST,
        windowmasker, taxonomy, refseq release), downloading ITIS taxonomy tables, and creating BioSQL databases.
        This class currently REQUIRES an instance of
        ProjectManagement to be used with the proj_mana parameter.

        :param email: The email of the user for using during the FTP.
        :type email: str
        :param driver: The driver used for creating the BioSQL databases.
        :type driver: str
        :param project: The name of the project.
        :type project: str or None
        :param project_path: A path used for standalone/basic project configuration.
        :type project_path: str or Path or None
        :param proj_mana: A configuration variable for connecting projects.
        :type proj_mana: ProjectManagement or None
        :param blast: Flag for BLAST-related database operations.
        :type blast: bool
        :param ftp_flag: A flag used if FTP connection is available or not.
        :type ftp_flag: bool
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

    @staticmethod
    def _prepare_child_strategies(
        child_strategies: tuple[dict[str, Any], ...],
        configure_flag: bool | None,
        archive_flag: bool | None,
        delete_flag: bool | None,
    ) -> tuple[dict[str, Any], ...]:
        """Copy inherited flags into child strategies without mutating inputs."""
        prepared_strategies = []
        for child_strategy in child_strategies:
            prepared_strategy = dict(child_strategy)
            if configure_flag:
                prepared_strategy["configure_flag"] = configure_flag
            if archive_flag:
                # The parent archive replaces child archive and delete actions.
                prepared_strategy["archive_flag"] = None
                prepared_strategy["delete_flag"] = None
            else:
                prepared_strategy["delete_flag"] = delete_flag
            prepared_strategies.append(prepared_strategy)
        return tuple(prepared_strategies)

    def _resolve_database_paths(
        self,
        database_path: str | Path | None,
        archive_path: str | Path | None,
    ) -> tuple[str | Path, str | Path]:
        """Resolve storage defaults while preserving explicit path objects."""
        resolved_database_path = database_path or str(self.user_db)
        resolved_archive_path = archive_path or str(self.user_archive)
        return resolved_database_path, resolved_archive_path

    def _append_archive_action(
        self,
        dispatcher: dict[str, Any],
        configuration: dict[str, Any],
        strategy_name: str,
        database_path: str | Path | None,
        archive_path: str | Path,
        delete_flag: bool | None,
        archive_option: str | None = None,
    ) -> None:
        """Append one archive action and its positionally matched configuration."""
        dispatcher.setdefault(strategy_name, []).append(self.db_mana_utils.archive)
        configuration.setdefault(strategy_name, []).append(
            {
                "database_path": database_path,
                "archive_path": archive_path,
                "option": archive_option or strategy_name,
                "delete_flag": delete_flag,
            }
        )

    def _build_leaf_strategy(
        self,
        *,
        strategy_name: str,
        configure_flag: bool | None,
        archive_flag: bool | None,
        delete_flag: bool | None,
        database_path: str | Path | None,
        archive_path: str | Path | None,
        configure_action: Callable[..., Any],
        configure_kwargs_factory: Callable[[], dict[str, Any]],
        archive_option: str | None = None,
        use_default_database_path: bool = True,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
        """Build one ordered leaf strategy from archive and configure actions."""
        dispatcher = OrderedDict({strategy_name: []})
        configuration = OrderedDict({strategy_name: []})
        resolved_database_path = database_path
        resolved_archive_path = archive_path or str(self.user_archive)
        if use_default_database_path and not resolved_database_path:
            resolved_database_path = str(self.user_db)

        if archive_flag:
            self._append_archive_action(
                dispatcher,
                configuration,
                strategy_name,
                resolved_database_path,
                resolved_archive_path,
                delete_flag,
                archive_option=archive_option,
            )
        if configure_flag:
            dispatcher[strategy_name].append(configure_action)
            configuration[strategy_name].append(configure_kwargs_factory())
        return dispatcher, configuration

    @staticmethod
    def _merge_strategy_results(
        dispatcher: dict[str, Any],
        configuration: dict[str, Any],
        child_results: tuple[tuple[dict[str, Any], dict[str, Any]], ...],
    ) -> None:
        """Merge paired child results in dispatch order."""
        for child_dispatcher, child_configuration in child_results:
            dispatcher.update(child_dispatcher)
            configuration.update(child_configuration)

    def get_strategy_dispatcher(
        self,
        db_config_strategy: dict[str, dict[str, Any]],
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
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
        strategy_methods = {
            "Full": self.full,
            "NCBI": self.NCBI,
            "NCBI_blast": self.NCBI_blast,
            "NCBI_blast_db": self.ncbi_blast_db,
            "NCBI_blast_windowmaskerfiles": self.ncbi_blast_windowmasker_files,
            "NCBI_pub_taxonomy": self.NCBI_pub_taxonomy,
            "NCBI_refseq_release": self.NCBI_refseq_release,
            "ITIS": self.itis,
            "ITIS_taxonomy": self.itis_taxonomy,
        }
        for strategy, strategy_kwargs in db_config_strategy.items():
            strategy_method = strategy_methods.get(strategy)
            if strategy_method is None:
                continue

            child_dispatcher, child_config = strategy_method(**strategy_kwargs)
            if strategy == "Full":
                strategy_dispatcher.clear()
                strategy_config.clear()
            strategy_dispatcher.update(child_dispatcher)
            strategy_config.update(child_config)

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

    def NCBI(
        self,
        NCBI_blast: dict[str, Any],
        NCBI_pub_taxonomy: dict[str, Any],
        NCBI_refseq_release: dict[str, Any],
        configure_flag: bool | None = None,
        archive_flag: bool | None = None,
        delete_flag: bool | None = None,
        database_path: str | Path | None = None,
        archive_path: str | Path | None = None,
        _path: str | Path | None = None,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
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
        database_path, archive_path = self._resolve_database_paths(
            database_path,
            archive_path,
        )
        (
            NCBI_blast,
            NCBI_pub_taxonomy,
            NCBI_refseq_release,
        ) = self._prepare_child_strategies(
            (NCBI_blast, NCBI_pub_taxonomy, NCBI_refseq_release),
            configure_flag,
            archive_flag,
            delete_flag,
        )

        if archive_flag:
            self._append_archive_action(
                ncbi_dispatcher,
                ncbi_config,
                "NCBI",
                database_path,
                archive_path,
                delete_flag,
            )

        nb_dispatcher, nb_config = self.NCBI_blast(**NCBI_blast)
        npt_dispatcher, npt_config = self.NCBI_pub_taxonomy(**NCBI_pub_taxonomy)
        nrr_dispatcher, nrr_config = self.NCBI_refseq_release(**NCBI_refseq_release)
        self._merge_strategy_results(
            ncbi_dispatcher,
            ncbi_config,
            (
                (nb_dispatcher, nb_config),
                (npt_dispatcher, npt_config),
                (nrr_dispatcher, nrr_config),
            ),
        )

        return ncbi_dispatcher, ncbi_config

    def NCBI_blast(
        self,
        NCBI_blast_db: dict[str, Any],
        NCBI_blast_windowmasker_files: dict[str, Any],
        configure_flag: bool | None = None,
        archive_flag: bool | None = None,
        delete_flag: bool | None = None,
        database_path: str | Path | None = None,
        archive_path: str | Path | None = None,
        _path: str | Path | None = None,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
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
        database_path, archive_path = self._resolve_database_paths(
            database_path,
            archive_path,
        )
        (
            NCBI_blast_db,
            NCBI_blast_windowmasker_files,
        ) = self._prepare_child_strategies(
            (NCBI_blast_db, NCBI_blast_windowmasker_files),
            configure_flag,
            archive_flag,
            delete_flag,
        )

        if archive_flag:
            self._append_archive_action(
                ncbi_blast_dispatcher,
                ncbi_blast_config,
                "NCBI_blast",
                database_path,
                archive_path,
                delete_flag,
            )

        nbd_dispatcher, nbd_config = self.ncbi_blast_db(**NCBI_blast_db)
        nbw_dispatcher, nbw_config = self.ncbi_blast_windowmasker_files(
            **NCBI_blast_windowmasker_files
        )
        self._merge_strategy_results(
            ncbi_blast_dispatcher,
            ncbi_blast_config,
            (
                (nbd_dispatcher, nbd_config),
                (nbw_dispatcher, nbw_config),
            ),
        )

        return ncbi_blast_dispatcher, ncbi_blast_config

    def ncbi_blast_db(
        self,
        configure_flag: bool | None = None,
        archive_flag: bool | None = None,
        delete_flag: bool | None = None,
        archive_path: str | Path | None = None,
        database_path: str | Path | None = None,
        _path: str | Path | None = None,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
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
        return self._build_leaf_strategy(
            strategy_name="NCBI_blast_db",
            configure_flag=configure_flag,
            archive_flag=archive_flag,
            delete_flag=delete_flag,
            database_path=database_path,
            archive_path=archive_path,
            configure_action=self.download_blast_database,
            configure_kwargs_factory=lambda: {"database_name": "refseq_rna"},
            archive_option="NCBI_blast",
            use_default_database_path=False,
        )

    def ncbi_blast_windowmasker_files(
        self,
        taxonomy_ids: Sequence[int],
        configure_flag: bool | None = None,
        archive_flag: bool | None = None,
        delete_flag: bool | None = None,
        archive_path: str | Path | None = None,
        database_path: str | Path | None = None,
        _path: str | Path | None = None,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
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
        # Preserve the established behavior of using project gene data here.
        _ = taxonomy_ids, _path
        return self._build_leaf_strategy(
            strategy_name="NCBI_blast_windowmasker_files",
            configure_flag=configure_flag,
            archive_flag=archive_flag,
            delete_flag=delete_flag,
            database_path=database_path,
            archive_path=archive_path,
            configure_action=self.download_windowmasker_files,
            configure_kwargs_factory=lambda: {
                "taxonomy_ids": self.gene_data.taxon_ids
            },
        )

    def NCBI_pub_taxonomy(
        self,
        configure_flag: bool | None = None,
        archive_flag: bool | None = None,
        delete_flag: bool | None = None,
        archive_path: str | Path | None = None,
        database_path: str | Path | None = None,
        _path: str | Path | None = None,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
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
        _ = _path
        return self._build_leaf_strategy(
            strategy_name="NCBI_pub_taxonomy",
            configure_flag=configure_flag,
            archive_flag=archive_flag,
            delete_flag=delete_flag,
            database_path=database_path,
            archive_path=archive_path,
            configure_action=self.download_ncbi_taxonomy_dump_files,
            configure_kwargs_factory=dict,
        )

    def NCBI_refseq_release(
        self,
        configure_flag: bool | None = None,
        archive_flag: bool | None = None,
        delete_flag: bool | None = None,
        upload_flag: bool | None = None,
        archive_path: str | Path | None = None,
        database_path: str | Path | None = None,
        collection_subset: str | None = None,
        seqtype: str | None = None,
        seqformat: str | None = None,
        file_list: list[str] | None = None,
        upload_number: int = 8,
        _path: str | Path | None = None,
        activate: str | Path | None = None,
        template_flag: bool | None = None,
        download_flag: bool | None = None,
        pbs_dict: dict[str, Any] | None = None,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
        """Build RefSeq download, taxonomy, and archive actions.

        The scheduler-backed BioSQL upload workflow has been retired. Its
        legacy arguments remain temporarily so existing configurations receive
        an explicit deprecation error instead of failing during argument
        parsing.

        :param configure_flag:  A flag that downloads refseq release files from
            NCBI.
        :type configure_flag:  bool.
        :param archive_flag:  A flag that archives refseq release files from NCBI.
        :type archive_flag:  bool.
        :param delete_flag:  A flag that deletes refseq release files from NCBI.
        :type delete_flag:  bool.
        :param upload_flag:  A retired flag that now raises a deprecation error.
        :type upload_flag:  bool.
        :param database_path:  User-supplied relative path to the databases.
        :type database_path:   str.
        :param archive_path:  User-supplied relative path to the archived
            databases.
        :type archive_path:   str.
        :param collection_subset: The collection subset of interest.
        :type collection_subset: str.
        :param seqtype: The type of sequence (rna, protein, genomic).
        :type seqtype: str.
        :param seqformat: The sequence format, usually ``gbff`` for GenBank Flat
            File.
        :type seqformat: str.
        :param file_list:  A deprecated upload argument retained for
            compatibility.
        :type file_list:  list.
        :param upload_number:  A deprecated upload argument retained for
            compatibility.
        :type upload_number:  int.
        :param activate:  A deprecated upload argument retained for compatibility.
        :type activate:  str.
        :param template_flag:  A flag that loads a BioSQL database with NCBI
            taxonomy data.
        :type template_flag:  bool.
        :param download_flag:  A flag that downloads RefSeq release files from
            NCBI's FTP site.
        :type download_flag:  bool.
        :param pbs_dict:  A deprecated upload argument retained for compatibility.
        :type pbs_dict:  dict.
        :return:  Paired action and configuration mappings for RefSeq data.
        :rtype:  tuple.
        """
        # These arguments remain in the signature only to give legacy callers
        # a precise error at the retired upload boundary.
        _ = file_list, upload_number, _path, activate, pbs_dict

        strategy_actions = OrderedDict(
            {"archive": [], "configure": [], "upload": []}
        )
        strategy_config = OrderedDict(
            {"archive": [], "configure": [], "upload": []}
        )
        nrr_dispatcher = OrderedDict(
            {"NCBI_refseq_release": strategy_actions}
        )
        nrr_config = OrderedDict({"NCBI_refseq_release": strategy_config})
        if upload_flag:
            raise OrthoEvolDeprecationWarning(
                "RefSeq BioSQL uploads depended on the retired SGE subsystem "
                "and are no longer supported."
            )

        dl_path = Path(self.database_path) / "NCBI" / "pub" / "taxonomy"
        dmp_file = dl_path / "nodes.dmp"
        if not archive_path:
            archive_path = str(self.user_archive)
        if not database_path:
            database_path = str(self.user_db)
        if archive_flag:
            nrr_dispatcher["NCBI_refseq_release"]["archive"].append(
                self.db_mana_utils.archive
            )
            nrr_config["NCBI_refseq_release"]["archive"].append(
                {
                    "database_path": database_path,
                    "archive_path": archive_path,
                    "option": "NCBI_refseq_release",
                    "delete_flag": delete_flag,
                }
            )
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

        return nrr_dispatcher, nrr_config

    def itis(
        self,
        ITIS_taxonomy: dict[str, Any],
        configure_flag: bool | None = None,
        archive_flag: bool | None = None,
        delete_flag: bool | None = None,
        database_path: str | Path | None = None,
        archive_path: str | Path | None = None,
        _path: str | Path | None = None,
    ) -> tuple[OrderedDict[str, Any], OrderedDict[str, Any]]:
        itis_dispatcher = OrderedDict()
        itis_config = OrderedDict()
        database_path, archive_path = self._resolve_database_paths(
            database_path,
            archive_path,
        )
        (ITIS_taxonomy,) = self._prepare_child_strategies(
            (ITIS_taxonomy,),
            configure_flag,
            archive_flag,
            delete_flag,
        )

        if archive_flag:
            self._append_archive_action(
                itis_dispatcher,
                itis_config,
                "ITIS",
                database_path,
                archive_path,
                delete_flag,
            )

        it_dispatcher, it_config = self.itis_taxonomy(**ITIS_taxonomy)
        self._merge_strategy_results(
            itis_dispatcher,
            itis_config,
            ((it_dispatcher, it_config),),
        )

        return itis_dispatcher, itis_config

    def itis_taxonomy(self, configure_flag=None, archive_flag=None, delete_flag=None, **kwargs):
        return {}, {}
