mpi module
----------

The mpi module will allow us to make use of the MCSR's PBS script
functionality via the ``mpi4py`` package. ## Description

Originally this was used for the FTP downloads. It is currently under
development for project wide use. The following need to be added: - [ ]
Start function - [ ] Function to split things into a nested list - [ ]
Function to create JSON file for each process - [ ] Remove functions
from the ftp2db script and put here - [ ] Add parameters for this module
- [ ] MP script - [ ] Script type (language) - [ ] Add a script for
calling these types (e.g. python3 script.py; R app.R) - [ ] Number of
processes - [ ] Nested list to be split - [ ] Flag for keeping default
PBS(.sh) script or generating a custom one - [ ] Look into the other
multiprocessing thing I linked on SLACK

Usage
-----

Usage will be more concisely described after this module is update.

A PBS script is called like this:

.. code:: bash

    $ qsub UPLOAD.sh

Here is what the current PBS script looks like:

.. code:: pbs

    #PBS -S /bin/bash
    #PBS -m bea
    #PBS -M rgilmore@umc.edu
    #PBS -l select=8:ncpus=1:mem=16gb -l place=free
    #PBS -l cput=24:00:00
    #PBS -l walltime=32:00:00
    #PBS -N robupload
    #PBS -o /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/Lib/Log/robupload.o${PBS_JOBID}
    #PBS -e /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/Lib/Log/robupload.e${PBS_JOBID}
    #PBS -j oe
    cd ${PBS_O_WORKDIR}
    rm /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/1_Databases/NCBI_Data/refseq/release/vertebrate_mammalian/robupload.o*
    rm /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/1_Databases/NCBI_Data/refseq/release/vertebrate_mammalian/robupload.e*
    mpiexec python /work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/1_Databases/NCBI_Data/refseq/release/vertebrate_mammalian/multi_dbupload.py
    echo "end"

This module is strictly a python driven module. ``mpi4py`` is used like
so:

.. code:: python

    from mpi4py import MPI

    # Get child process information
    comm = MPI.COMM_WORLD

    # The rank is unique to each process.
    # For 8 parallel processes there will be a rank 0-7
    rank = comm.Get_rank()
    # Currently unused variable in my scripts
    size = comm.Get_size()
    machine = platform.node()

If another type of programming language is needed then ## Tests

??
