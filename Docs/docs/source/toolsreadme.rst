Tools Documentation
===================

These tools were created by or modified our team to aid with the
Orthologs package.

We've incorporated tools for bash with pybasher, qsub tools for use with
pbs, a pandoc script for converting docx files to markdown formats,
multiprocessing in multiprocess, and a ftp module that aids in
downloading files from NCBI's ftp repository.

Examples
--------

Download NCBI databases with our NCBI FTP Client
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from Datasnakes.Tools.ftp import NcbiFTPClient

    ncbiftp = NcbiFTPClient(email='somebody@gmail.com')
    ncbiftp.getblastdb(database_name='refseq_rna')

List all subdirectories in a NCBI FTP Path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python


    ncbiftp.listdirectories(path='/blast/db/')
    Out[54]: ['FASTA', 'cloud']

Utilize multiprocessing to speed up your code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from Datasnakes.Tools import Multiprocess


    def printwords(word):
        print(word)


    words = ['bae', 'luh', 'cuh']

    if __name__ == '__main__':
        mp = Multiprocess()
        mp.map2function(printwords, words)

Integrate logging in a simple and quick way
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from Datasnakes.Tools import LogIt

    # Set up your loggers
    logit = LogIt()

    # Log to one file
    logfile = 'test.log'

    test1 = logit.default('test1 log', logfile)

    # Start logging
    test1.info('hi')

    # Shutdown logging without deleting the logfile
    logit.shutdown()

Additional Documentation
------------------------

Check the specific modules for more detailed readmes and examples of
using the tools with this package.
