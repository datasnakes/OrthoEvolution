## biosql directory

This is the biosql package for our project.  It will automate upload of the
sqlite3 biosql schema and the NCBI taxonomy data to a template database
once the package is installed via pip.
> Much of this README is referenced from [biopython's](http://biopython.org)
and [biosql's](http://biosql.org) documentation site.

## Description

#### BioSQL

BioSQL is a generic unifying schema for storing sequences from different sources,
for instance Genbank or Swissprot.

BioSQL is meant to be a common data storage layer supported by all the
different Bio* projects, Bioperl, Biojava, Biopython, and Bioruby.
Entries stored through an application written in, say, Bioperl could
be retrieved by another written in Biojava.

There are currently four different RDBMSs supported: MySQL,
PostgreSQL, Oracle, and most recently SQLite. The MySQL schema
DDL is in sql/biosqldb-mysql.sql, the PostgreSQL schema is in
sql/biosqldb-pg.sql, and the Oracle schema is in multiple files
in the directory sql/biosql-ora, and the SQLite schema is in
sql/biosqldb-sqlite.sql. In order to instantiate the schema, feed
the respective file or files to your SQL shell (e.g. mysql for
MySQL, and psql for PostgreSQL).

>Note the SQLite schema is new, and is not yet supported by all the
>Bio* projects (at the time of writing, just Biopython supports this).

#### PhyloDB

PhyloDB is an extension of the BioSQL schema.  It lies on top of an
existing BioSQL database and gives the added functionality of storing
and querying data from common phylogenetic file types (e.g. NEXUS and Newick
files).
>Currently only available for MySQL and PGSQL.

## Usage

#### Custom dependencies

- Python3.6
- biopython
- sqlite3

    ```bash
    $ wget http://www.sqlite.org/sqlite-autoconf-3070603.tar.gz
    $ tar xvfz sqlite-autoconf-3070603.tar.gz
    $ cd sqlite-autoconf-3070603
    $ ./configure
    $ make
    $ make install
    ```
- perl (NCBI taxonomy loading)
    ```bash
    $ sudo apt-get install perl
    ```
- MySQL (future for PhyloDB)

#### Schema Loading

The BioSQL schema is distributed separately from the language
bindings. Use git clone to get a local copy:

```bash
$ git clone https://github.com/biosql/biosql.git
```
Create a database, which we'll call 'biosql', in the data instance:

* For SQLITE3 do:
    ```bash
    $ sqlite3 biosql.db
    ```
* For MySQL do:
    ```bash
    $ mysqladmin -u root create biosql
    ```

To load the schema, use the appropiate SQL dialect in
biosql-schema/sql.

* For SQLITE3 do:
    ```bash
    $ sqlite3 biosqldb < biosqldb-sqlite.sql
    ```
* For mysql do:
    ```bash
    $ mysql -u root biosql < biosqldb-mysql.sql
    ```

#### NCBI taxonomy loading

Loading NCBI taxonomy requires perl and a BioSQL database loaded with the
BioSQL schema.

The BioSQL package includes a perl script under scripts/load_ncbi_taxonomy.pl
to download and update the taxonomy tables. The script should be able to download
the files it needs from the NCBI taxonomy FTP site automatically.

To update the NCBI taxonomy, change to the scripts subdirectory from the
BioSQL repository, then copy and paste the following command:

```bash
$ ./load_ncbi_taxonomy.pl --dbname biosql.db --driver sqlite3 --dbuser root --download true
```

#### Loading sequence data

For our uses loading in sequence data requires biopython.
BioSQL lets us define named “sub” databases or “namespaces” within the
single SQL database (which we called biosql.db earlier). For this example,
lets create a one for some orchid sequences:
```python
from BioSQL import BioSeqDatabase
server = BioSeqDatabase.open_database(driver="sqlite3", db="biosql.db")
db = server.new_database("orchids", description="Just for testing")
...
```
The *commit* call tells the database to save the changes so far (commit
the SQL transaction). It is up to you to decide when to commit the SQL
transaction(s), and/or rollback changes, rather than having Biopython
try and decide for you and risk getting it wrong. See *Explicit is
better than implicit* ([The Zen of
Python](http://www.python.org/dev/peps/pep-0020/)).

There should now be a single row in the *biodatabase* table for our new
orchid namespace.

Now lets continue with this example by loading some data from NCBI using
biopython's Entrez module.

```
...
from Bio import Entrez
from Bio import SeqIO
handle = Entrez.efetch(db="nuccore", id="6273291,6273290,6273289", rettype="gb", retmode="text")
count = db.load(SeqIO.parse(handle, "genbank"))
print "Loaded %i records" % count
server.commit()
```
Again, you must explicitly call commit to record the SQL transaction which is otherwise left pending.

The db.load() function should have returned the number of records loaded
(three in this example), and again have a look in the database and you should see new rows in several tables.

#### Extracting sequence data

This continues from the previous example, where we loaded three records into an orchids database (namespace):
```python
from BioSQL import BioSeqDatabase
server = BioSeqDatabase.open_database(driver="sqlite3", db="biosql.db")
db = server["orchids"]
for identifier in ['6273291', '6273290', '6273289'] :
    seq_record = db.lookup(gi=identifier)
    print seq_record.id, seq_record.description[:50] + "..."
    print "Sequence length %i," % len(seq_record.seq)
```

#### Deleting a subdatabase
As mentioned above, BioSQL lets us define named "sub" databases (aka namespaces)
within the single SQL database (which we called bioseqdb). In the previous example,
we created a sub-database for some orchid sequences. The following code will delete
the orchid database (and all the records in it):

```python
from BioSQL import BioSeqDatabase
server = BioSeqDatabase.open_database(driver="sqlite3", db="biosql.db")
server.remove_database("orchids")
server.commit()
```