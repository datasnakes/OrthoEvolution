import pkg_resources
from pathlib import Path
import os
import sys
import subprocess
import shutil
from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Datasnakes.Tools.logit import LogIt
from Datasnakes.Manager.BioSQL.biosql_repo import sql
from Datasnakes.Manager.BioSQL.biosql_repo import scripts as sql_scripts
from Datasnakes.Manager.management import ProjectManagement
from Datasnakes.Orthologs.utils import attribute_config
from Datasnakes.Tools.streamieo import StreamIEO


class BaseBioSQL(object):
    # TODO-ROB:  Organize the BioSQL files by driver/RDBMS
    # TODO-ROB:  Add functionality for database_type="biosqldb"
    def __init__(self, database_name, template_name="", project=None, project_path=None, proj_mana=ProjectManagement, **kwargs):
        """
        This is the base BioSQL class.  It provides a framework for uploading schemas, loading taxonomy data from NCBI
        and ITIS using the BioSQL perl scripts and .sql schema files provided by the BioPython package.  We have created
        a modified version of the BioSQL scripts in our package, which can be found on GitHub.  Taxonomy data can be
        found at:
            NCBI:  ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy
            ITIS:  http://www.itis.gov/downloads/

        :param database_name:  The name of the database.
        """
        # Logging setup
        self.biosqllog = LogIt().default(logname="BioSQL", logfile=None)
        self.biosqlstream = StreamIEO(logname="BioSQL")

        # Load relative and absolute paths to scripts in the BioSQL module
        self.scripts = pkg_resources.resource_filename(sql_scripts.__name__, "")
        self.ncbi_taxon_script = pkg_resources.resource_filename(sql_scripts.__name__, "load_ncbi_taxonomy.pl")
        self.itis_taxon_script = pkg_resources.resource_filename(sql_scripts.__name__, "load_itis_taxonomy.pl")
        self.database_name = database_name

        # Configuration of class attributes for Project Management.
        if project_path and project:
            self.project_path = Path(project_path) / Path(project)

        if proj_mana:
            add_self = attribute_config(self, composer=proj_mana, checker=ProjectManagement, project=project, project_path=project_path)
            for var, attr in add_self.__dict__.items():
                setattr(self, var, attr)
            self.template_rel_path = self.user_index
            self.template_abs_path = self.template_rel_path / Path(template_name)
            self.database_rel_path = self.user_db
            self.database_abs_path = self.database_rel_path / Path(self.database_name)
        else:
            self.project_path = Path(project_path) / Path(project)
            self.template_rel_path = self.project_path / Path('index')
            self.template_abs_path = self.template_rel_path / Path(template_name)
            self.database_rel_path = Path(project_path) / Path('databases')
            self.database_abs_path = self.database_rel_path / Path(self.database_name)

    def configure_new_database(self, cmd, schema_file=None):
        """
        This script is a framework for loading the various schemas, the NCBI taxonomy (biosql-db), and the ITIS
        taxonomy (phylo-db) into a database.

        :param cmd:  The bash command to use.
        :param schema_file:  The schema file for creating a BioSQL or PhyloDB
        :return:  Returns the Output and the Error messages.
        """
        # Build the command
        # To upload a schema, a schema file will be necessary.
        if schema_file:
            _ = "Schema"
            cmd = "%s < %s" % (cmd, schema_file)
        # But not for uploading taxonomy data (NCBI/ITIS)
        else:
            _ = "Taxonomy"
            cmd = cmd

        # Run the perl scripts or sql scripts
        self.biosqlstream.streamer(cmd)

    def create_executable_scripts(self):
        """
        Changes the permissions of the BioSQL perl scripts, so that they are executable from the command line.
        """
        # Set up the permissions for the BioSQL Perl scripts
        biosql_scripts = self.scripts
        for file in os.listdir(biosql_scripts):
            print(file)
            if '.pl' in str(file):
                script_path = os.path.join(biosql_scripts, file)
                print(script_path)
                os.chmod(script_path, mode=0o755)


class SQLiteBioSQL(BaseBioSQL):
    def __init__(self, database_name=None, proj_mana=ProjectManagement, template_name="Template-BioSQL-SQLite.db", upload_path=None, upload_list=None, **kwargs):
        """
        This class inherits the BaseBioSQL class.  It uses the base methods to load schema, load taxonomy (NCBI),
        and create/copy template SQLite databases loaded with biosql schema and/or taxonomy data.

        :param database_name:  The name of the database.
        """
        super().__init__(template_name=template_name, database_name=database_name, proj_mana=proj_mana, **kwargs)
        self.driver = "SQLite"
        self.schema_cmd = "sqlite3 %s -echo"
        self.schema_file = "biosqldb-sqlite.sql"
        self.taxon_cmd = "%s --dbname %s --driver %s --download false --directory %s"

        self.upload_path = upload_path
        self.upload_list = upload_list

    def load_sqlite_schema(self):
        """
        Loads an SQLite biosql schema into a database file.
        """
        # Build the command
        schema_file = pkg_resources.resource_filename(sql.__name__, self.schema_file)
        schema_cmd = self.schema_cmd % str(self.template_abs_path)
        # Run the bash command
        self.configure_new_database(schema_cmd, schema_file)

        # TODO-ROB: Parse output and error
        # TODO-ROB:  Make sure the .db file doesn't already exist

    def load_sqlite_taxonomy(self):
        """
        Load an SQLite biosql database with taxonomy information.  This will only work for NCBI.  There is not any
        support for the SQLite PhyloDB.
        """
        # Build the command
        ncbi_taxon_dump_path = self.database_rel_path / Path("NCBI") / Path('pub') / Path('taxonomy')
        taxon_cmd = self.taxon_cmd % (self.ncbi_taxon_script, str(self.template_abs_path), self.driver, str(ncbi_taxon_dump_path))
        # Run the bash command
        self.configure_new_database(taxon_cmd)
        # TODO-ROB: Parse output and error
        # TODO-ROB:  Make sure the .db file doesn't already exist

    def create_template_database(self):
        """
        Creates a template database by uploading SQLite schema and NCBI taxonomy.
        :return:
        """
        # Create a template if it doesn't exits.
        if not self.template_abs_path.is_file():
            self.load_sqlite_schema()
            self.create_executable_scripts()
            # TODO-ROB:  Download the file manually and then use qsub job
            self.load_sqlite_taxonomy()
        else:
            self.biosqllog.warning("The template, %s, already exists." % self.template_abs_path)

    def copy_template_database(self, destination):
        """
        This method copies a template sqlite biosql database.

        :param destination:  The path to copy the template into
        :return:  A new copy of the biosql database.
        """
        # If the template doesn't exists, then create it.
        if not self.template_abs_path.is_file():
            self.create_template_database()
        # Copy the template into a new folder.
        dest_abs_path = Path(destination)
        self.biosqllog.warn('Copying Template BioSQL Database...  This may take a few minutes...')
        shutil.copy2(str(self.template_abs_path), str(dest_abs_path))

    def upload_files(self, seqtype, filetype, new_db=False):
        db_name = Path(self.database_abs_path.stem + '_' + seqtype + self.database_abs_path.suffix)
        db_abs_path = self.database_rel_path / db_name

        # Make sure a BioSQL-SQLite database exists
        if self.database_abs_path.is_file():
            if not db_abs_path.is_file():
                self.database_abs_path.rename(target=db_abs_path)
            pass
        elif new_db:
            self.copy_template_database(dest_path=self.database_rel_path, dest_name=db_name)
        else:
            raise FileNotFoundError("Database not found: %s\mPlease create a BioSQL-SQLite database." % self.database_abs_path)

        # Parse the upload list and upload the files to the BioSQL-SQLite database.
        for file in self.upload_list:
            abs_upload_path = Path(str(self.upload_path)) / Path(file)

            # Make a connection with the BioSQL database
            try:
                server = BioSeqDatabase.open_database(driver=self.driver.lower(), db=str(db_abs_path))
                self.biosqllog.info("Server Connected.")
                pass
            except:
                self.biosqllog.warn("The Server did not Connect.  Check the to make sure %s exists." % self.database_abs_path)
                raise FileNotFoundError

            # See if the sub database exists (rna, protein, or genomic)
            try:
                if seqtype not in server.keys():
                    server.new_database(seqtype)
                    self.biosqllog.info("New Sub-Database created, %s, for %s." % (seqtype, db_abs_path))
                # Connect to the sub database
                sub_db = server[seqtype]

                count = sub_db.load(SeqIO.parse(abs_upload_path, filetype))
                self.biosqllog.info("%s loaded with %s %s files" % (db_name, count, filetype))
                server.commit()
                self.biosqllog.warn("Server committed.")
                t_count = t_count + count
                self.biosqllog.info("The server has not loaded a total of %s files." % t_count)
                # TODO-ROB:  Add something to do with time here.
            except:
                self.biosqllog.critical("Unable to load the database...")
                server.rollback()
                try:
                    del server[sub_db]
                    self.biosqllog.critical("%s sub database deleted from %s.  All of the info will be lost." % (sub_db, db_abs_path))
                    server.commit()
                    self.biosqllog.critical("Server committed")
                except:
                    raise
                raise


class MySQLBioSQL(BaseBioSQL):

    def __init__(self, database_path, database_name="Template-BioSQL-MySQL.db"):
        """
        This class inherits the BaseBioSQL class.  It uses the base methods to load schema, load taxonomy (NCBI, ITIS),
        and create/copy template MySQL databases loaded with biosql schema and/or taxonomy data.  The MySQL driver
        for BioPython's BioSQL databases is the most developed.  It can utilize the PhyloDB schema extension and scripts
        as well as the standard BioSQL schema and scripts.

        :param database_path:  The relative path to the database.
        :param database_name:  The name of the database.
        """
        super().__init__(database_name=database_name, database_path=database_path, driver="MySQL")
        self.schema_cmd = "sqlite3 %s -echo"
        self.schema_file = "biosqldb-mysql.sql"
        self.taxon_cmd = "%s --dbname %s --driver %s --download true"
