import pkg_resources
from pathlib import Path
import os
import subprocess
import shutil
from Datasnakes.Tools.logit import LogIt
from Datasnakes.Manager.BioSQL.biosql_repo import sql
from Datasnakes.Manager.BioSQL.biosql_repo import scripts as sql_scripts


class BaseBioSQL(object):
    # TODO-ROB:  Organize the BioSQL files by driver/RDBMS
    # TODO-ROB:  Add functionality for database_type="biosqldb"
    def __init__(self, database_name, database_path, driver):
        """
        This is the base BioSQL class.  It provides a framework for uploading schemas, loading taxonomy data from NCBI
        and ITIS using the BioSQL perl scripts and .sql schema files.  We have created a modified version of the
        BioSQL file system in our package, which can be found on GitHub.
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
    def __init__(self, database_name="Template-BioSQL-SQLite.db"):
        """
        This class inherits the BaseBioSQL class.  It uses the base methods to load schema, load taxonomy (NCBI, ITIS),
        and create/copy template sqlite databases loaded with biosql schema and/or taxonomy data.
        :param database_name:  The name of the database.
        :param template:  The standard template name.
        """
        super().__init__(database_name=database_name, database_path="", driver="SQLite")
        self.schema_cmd = "sqlite3 %s -echo"
        self.schema_file = "biosqldb-sqlite.sql"
        self.taxon_cmd = "%s --dbname %s --driver %s --download true"

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
            self.biosqllog.warning("The template, %s, already exists." % db_path)

    def copy_template_database(self, dest_path, dest_name=""):
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


class MySQLBioSQL(BaseBioSQL):

    def __init__(self):
        super().__init__(database_name="", database_path="", driver="")
        pass
