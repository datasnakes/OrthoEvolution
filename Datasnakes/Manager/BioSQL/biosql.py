import pkg_resources
from pathlib import Path
import os
import subprocess
from Datasnakes.Tools.logit import LogIt
from Datasnakes.Manager.BioSQL.biosql import sql
from Datasnakes.Manager.BioSQL.biosql import scripts as sql_scripts


class BaseBioSQL(object):
 # driver "mysql", "Pg", "Oracle", "SQLite"
    def __init__(self, database_name, database_type, driver):
        self.database_name = database_name
        self.database_type = database_type
        self.driver = driver
        self.biosqllog = LogIt().default(logname="BioSQL", logfile=None)
        pass

    def load_biosql_schema(self, cmd, schema_file):
        biosql_cmd = "%s < %s" % (cmd, schema_file)
        biosql_load = subprocess.Popen([biosql_cmd], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True,
                                       encoding='utf-8')
        error = biosql_load.stderr.readlines()
        out = biosql_load.stdout.readlines()

        self.biosqllog.info("Error: " + str(error))
        self.biosqllog.info("Out: " + str(out))

    def load_ncbi_taxonomy(self):
        pass


class SQLiteBioSQL(BaseBioSQL):
    def __init__(self, database_name):
        super().__init__(database_name=database_name, driver="SQLite")
        self.cmd = "sqlite3 %s -echo" % database_name
        self.schema_file = "biosqldb-sqlite.sql"

    def sqlite_schema(self):
        self.load_biosql_schema(self.cmd, self.schema_file)
