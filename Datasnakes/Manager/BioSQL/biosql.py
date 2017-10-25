import pkg_resources
from pathlib import Path
import os
import subprocess
from Datasnakes.Tools.logit import LogIt
from Datasnakes.Manager.BioSQL.biosql_repo import sql
from Datasnakes.Manager.BioSQL.biosql_repo import scripts as sql_scripts


class BaseBioSQL(object):
    # TODO-ROB:  Organize the BioSQL files by driver/RDBMS
 # driver "mysql", "Pg", "Oracle", "SQLite"
    def __init__(self, database_name, database_type, driver):
        self.database_name = database_name
        self.database_type = database_type
        self.driver = driver
        self.biosqllog = LogIt().default(logname="BioSQL", logfile=None)

        self.scripts = pkg_resources.resource_filename(sql_scripts.__name__, "")
        self.ncbi_taxon_script = pkg_resources.resource_filename(sql_scripts.__name__, "load_ncbi_taxonomy.pl")
        pass

    def load_biosql_schema(self, cmd, schema_file):
        schema_cmd = "%s < %s" % (cmd, schema_file)
        schema_load = subprocess.Popen([schema_cmd], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True,
                                       encoding='utf-8')
        error = schema_load.stderr.readlines()
        out = schema_load.stdout.readlines()

        self.biosqllog.info("Schema-Error: " + str(error))
        self.biosqllog.info("Schema-Out: " + str(out))
        return error, out

    def load_ncbi_taxonomy(self, cmd):
        # ./load_ncbi_taxonomy.pl --dbname bioseqdb --driver mysql --dbuser root --download true
        taxon_cmd = cmd
        taxon_load = subprocess.Popen([taxon_cmd], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True,
                                      encoding='utf-8')

        error = taxon_load.stderr.readlines()
        out = taxon_load.stdout.readlines()

        self.biosqllog.info("Taxon-Error: " + str(error))
        self.biosqllog.info("Taxon-Out: " + str(out))
        return error, out
        pass

    def create_executable_scripts(self):
        # Set up the permissions for the BioSQL Perl scripts
        biosql_scripts = self.scripts
        for file in os.listdir(biosql_scripts):
            print(file)
            if '.pl' in str(file):
                script_path = os.path.join(biosql_scripts, file)
                print(script_path)
                os.chmod(script_path, mode=0o755)


class SQLiteBioSQL(BaseBioSQL):
    def __init__(self, database_name, database_type):
        super().__init__(database_name=database_name, database_type=database_type, driver="SQLite")
        self.schema_cmd = "sqlite3 %s -echo" % database_name
        self.schema_file = "biosqldb-sqlite.sql"

        self.taxon_cmd = "%s --dbname %s --driver %s --download true" % \
                         (self.ncbi_taxon_script, database_name, self.driver)

    def sqlite_schema(self):
        schema_file = pkg_resources.resource_filename(sql.__name__, self.schema_file)
        error, out = self.load_biosql_schema(self.schema_cmd, schema_file)
        # TODO-ROB: Parse output and error

    def sqlite_taxonomy(self):
        error, out = self.load_ncbi_taxonomy(self.taxon_cmd)
        # TODO-ROB: Parse output and error