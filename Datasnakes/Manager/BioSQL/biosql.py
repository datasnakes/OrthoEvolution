import pkg_resources
from pathlib import Path
import os
import subprocess
import shutil
from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Datasnakes.Tools.logit import LogIt
from Datasnakes.Manager.BioSQL.biosql_repo import sql
from Datasnakes.Manager.BioSQL.biosql_repo import scripts as sql_scripts


class BaseBioSQL(object):
    # TODO-ROB:  Organize the BioSQL files by driver/RDBMS
    # TODO-ROB:  Add functionality for database_type="biosqldb"
    def __init__(self, database_name, database_path, driver):
        """
        This is the base BioSQL class.  It provides a framework for uploading schemas, loading taxonomy data from NCBI
        and ITIS using the BioSQL perl scripts and .sql schema files provided by the BioPython package.  We have created
        a modified version of the BioSQL scripts in our package, which can be found on GitHub.  Taxonomy data can be
        found at:
            NCBI:  ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy
            ITIS:  http://www.itis.gov/downloads/

        :param database_name:  The name of the database.
        :param driver:  The driver type.  "MySQL" (stable), "SQLite", "PostGRE"
        """
        # Logging setup
        self.biosqllog = LogIt().default(logname="BioSQL", logfile=None)

        # Parameter Attributes
        self.database_name = database_name
        self.database_path = database_path
        self.db_abs_path = Path(self.database_path) / Path(self.database_name)
        self.driver = driver

        # Load relative and absolute paths to scripts in the BioSQL module
        self.scripts = pkg_resources.resource_filename(sql_scripts.__name__, "")
        self.ncbi_taxon_script = pkg_resources.resource_filename(sql_scripts.__name__, "load_taxonomy.pl")
        self.itis_taxon_script = pkg_resources.resource_filename(sql_scripts.__name__, "load_itis_taxonomy.pl")

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
        _loader = subprocess.Popen([cmd], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, encoding='utf-8')
        out = _loader.stdout.readlines()
        error = _loader.stderr.readlines()

        self.biosqllog.info(_ + "-Error: " + str(error))
        self.biosqllog.info(_ + "-Out: " + str(out))
        return out, error

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
    def __init__(self, database_path=None, database_name="Template-BioSQL-SQLite.db", upload_path=None, upload_list=None):
        """
        This class inherits the BaseBioSQL class.  It uses the base methods to load schema, load taxonomy (NCBI),
        and create/copy template SQLite databases loaded with biosql schema and/or taxonomy data.

        :param database_path:  The relative path to the database.
        :param database_name:  The name of the database.
        """
        super().__init__(database_name=database_name, database_path=database_path, driver="SQLite")
        self.schema_cmd = "sqlite3 %s -echo"
        self.schema_file = "biosqldb-sqlite.sql"
        self.taxon_cmd = "%s --dbname %s --driver %s --download true"

        self.upload_path = upload_path
        self.upload_list = upload_list

    def load_sqlite_schema(self):
        """
        Loads an SQLite biosql schema into a database file.
        """
        # Build the command
        schema_file = pkg_resources.resource_filename(sql.__name__, self.schema_file)
        schema_cmd = self.schema_cmd % str(self.db_abs_path)
        # Run the bash command
        out, error = self.configure_new_database(schema_cmd, schema_file)
        return out, error
        # TODO-ROB: Parse output and error

    def load_sqlite_taxonomy(self):
        """
        Load an SQLite biosql database with taxonomy information.  This will only work for NCBI.  There is not any
        support for the SQLite PhyloDB.
        """
        # Build the command
        taxon_cmd = self.taxon_cmd % (self.ncbi_taxon_script, str(self.db_abs_path), self.driver)
        # Run the bash command
        out, error = self.configure_new_database(taxon_cmd)
        return out, error
        # TODO-ROB: Parse output and error

    def create_template_database(self):
        """
        Creates a template database by uploading SQLite schema and NCBI taxonomy.

        :param db_path:  The relative path of the database.
        :return:
        """
        # Create a template if it doesn't exits.
        if not self.db_abs_path.is_file():
            self.load_sqlite_schema()
            self.create_executable_scripts()
            # TODO-ROB:  Download the file manually and then use qsub job
            self.load_sqlite_taxonomy()
        else:
            self.biosqllog.warning("The template, %s, already exists." % self.db_abs_path)

    def copy_template_database(self, dest_path, dest_name):
        """
        This method copies a template sqlite biosql database.

        :param db_path:  The relative path of the template database.
        :param dest_path:  The path to copy the template into
        :param dest_name:  The name of the new database file.
        :return:  A new copy of the biosql database.
        """
        # If the template doesn't exists, then create it.
        if not self.db_abs_path.is_file():
            self.create_template_database()
        # Copy the template into a new folder.
        dest_abs_path = Path(dest_path) / Path(dest_name)
        self.biosqllog.warn('Copying Template BioSQL Database...  This may take a few minutes...')
        shutil.copy2(str(self.db_abs_path), str(dest_abs_path))

    def upload_files(self, seqtype, filetype, new_db=False):
        db_path = self.db_abs_path.parent
        db_name = Path(self.db_abs_path.stem + '_' + seqtype + self.db_abs_path.suffix)
        db_abs_path = db_path / db_name

        # Make sure a BioSQL-SQLite database exists
        if self.db_abs_path.is_file():
            if not db_abs_path.is_file():
                self.db_abs_path.rename(target=db_abs_path)
            pass
        elif new_db:
            self.copy_template_database(dest_path=db_path, dest_name=db_name)
        else:
            raise FileNotFoundError("Database not found: %s\mPlease create a BioSQL-SQLite database." % self.db_abs_path)

        # Parse the upload list and upload the files to the BioSQL-SQLite database.
        for file in self.upload_list:
            abs_upload_path = Path(str(self.upload_path)) / Path(file)

            # Make a connection with the BioSQL database
            try:
                server = BioSeqDatabase.open_database(driver=self.driver, db=str(db_abs_path))
                self.biosqllog.info("Server Connected.")
                pass
            except:
                self.biosqllog.warn("The Server did not Connect.  Check the to make sure %s exists." % self.db_abs_path)
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
