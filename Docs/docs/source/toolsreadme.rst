Tools Documentation
===================

The Tools module is a collection of often used classes or functions that
either enhance our other modules and create reusable functions to be
used in various modules.

We've incorporated tools for sge tools for use with pbs, a pandoc script
and class for converting docx files to markdown formats, multiprocessing
in multiprocess, and a ftp module that aids in downloading files from
NCBI's ftp repository.

Examples
--------

Take a look at the examples below to get an idea of how to incorporate
these tools in your project and how we use these tools in our project.

Download NCBI databases with our NCBI FTP Client
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools.ftp import NcbiFTPClient

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

    from OrthoEvol.Tools import Multiprocess


    def printwords(word):
        print(word)


    words = ['bae', 'luh', 'cuh']

    if __name__ == '__main__':
        mp = Multiprocess()
        mp.map2function(printwords, words)

Integrate logging in a simple and quick way
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools import LogIt

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
