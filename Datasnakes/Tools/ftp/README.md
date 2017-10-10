FTP (File Transfer Protocol) Documentation
=============================================
The `ftp` module is geared towards making it easier to interface with [NCBI's
FTP repository](ftp://ftp.ncbi.nlm.nih.gov).

More specifically, we provide a way to easily find and list directories and their
respective contents as well as to download blast databases and other databases
for use with the Orthologs package.

Examples
-----
These tools are optimized to be used together (very little work to do that), but can also be used singularly.


#### Blastdb Download Example

This is a simple example of using some of the modules.

``` python
from Datasnakes.Tools.ftp import NcbiFTPClient

ncbiftp = NcbiFTPClient(email='somebody@gmail.com')
ncbiftp.getblastdb(database_name='refseq_rna')

```

#### List all directories in a path
```python

ncbiftp.listdirectories('/blast/db/')
Out[54]: ['FASTA', 'cloud']
```

#### List all files in a path
```python

ncbiftp.listfiles('/blast/db/')
```

:exclamation: Notes
-------------------
Check the [NCBI README](NCBIREADME.md) for information about the preformatted blast databases that we use
and suggest you use. We also provide an easy way to download them which is referenced in the above example.


If you're using Linux and do not want to use threading to download ftp databases,
you can look at [this]() standalone script.