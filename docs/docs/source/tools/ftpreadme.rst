FTP (File Transfer Protocol) Documentation
==========================================

The ``ftp`` module is geared towards making it easier to interface with
`NCBI's FTP repository <ftp://ftp.ncbi.nlm.nih.gov>`__.

More specifically, we provide a way to easily find and list directories
and their respective contents as well as to download blast databases and
other databases for use with the Orthologs package. We have implemented
database downloading with threading which is the safest way to implement
this cross-platform.

We also provide a parallel module which can be used in conjunction with
the ``NcbiFTPClient`` to download files or databases much quicker if
your system can handle that.

If you're using Linux or a supercomputer and do not want to use
threading to download ftp databases, you can look at `this cli
script <https://github.com/datasnakes/OrthoEvolution/blob/master/Examples/standalone-scripts/ncbi-download.py>`__.

Examples
--------

Blastdb Download Example
^^^^^^^^^^^^^^^^^^^^^^^^

This is a simple example of using some of the modules.

.. code:: python

    from OrthoEvol.Tools import NcbiFTPClient

    ncbiftp = NcbiFTPClient(email='somebody@gmail.com')
    ncbiftp.getblastdb(database_name='refseq_rna')

Windowmasker files Download Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    from OrthoEvol.Tools import NcbiFTPClient
    import os

    ids = ['9544', '9606']

    ncbiftp = NcbiFTPClient(email='somebody@gmail.com')
    ncbiftp.getwindowmaskerfiles(taxonomy_ids=ids, download_path=os.getcwd())

Refseq Release Download Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    from OrthoEvol.Tools import NcbiFTPClient
    import os

    ncbiftp = NcbiFTPClient(email='somebody@gmail.com')
    ncbiftp.getrefseqrelease(taxon_group='vertebrate_mammalian', seqtype='rna', seqformat='gbff', download_path=os.getcwd())

List all directories in a path
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python


    ncbiftp.listdirectories(path='/blast/db/')
    Out[54]: ['FASTA', 'cloud']

List all files in a path
^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python


    ncbiftp.listfiles(path='/blast/db/')

List all files in the current working directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python


    # The default path is ftp.pwd() or the current directory
    ncbiftp.listfiles()

Notes
-----

Check the `NCBI README <NCBIREADME.md>`__ for information about the
preformatted blast databases that we use and suggest you use. We also
provide an easy way to download them which is referenced in the above
example.
