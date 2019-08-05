SGE Documentation
=================

Collection of tools for using PBS, a job scheduler for high-performance
computing environments on SGE. The command is usually ``qsub <options>``
on most systems.

Usage & Examples
----------------

The main class under sge is ``SGEJob``, which provides functionality to
use the job sceduling system on a high performance computing (HPC)
cluster.

The ``Qstat`` class is also available for parsing the output of the
``qstat`` command.

The class currently provides a template, ``temp.pbs``, file to be
modified and used when submitting a job as well as default job
attributes.

Using SGEJob with Multiprocess
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools.sge import SGEJob

Submitting multiple jobs
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools.sge import SGEJob

Get Job Info
~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools.sge import SGEJob

Running a simple job
~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools.sge import SGEJob

    myjob = SGEJob(email_address='shutchins2@umc.edu')

    code = "test.py"
    myjob.submit_pycode(code)

Software Dependencies
---------------------

Ensure that you have at least pbs version ``14.1.0``

Thanks
------

Thanks to [@jfeala](https://github.com/jfeala) for his work on Luigi's
SGEJobTask.
