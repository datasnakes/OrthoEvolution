sge Documentation
=================

Collection of tools for using PBS, a job scheduler for high-performance
computing environments on SGE. The command is usually ``qsub <options>``
on most systems.

This module also incorporates ``qstat <options>``.

Usage & Examples
----------------

The main class under sge is ``SGEJob``. Some functions are
``import_temp``, which allows the user to import a preformatted template
pbs script or python script and use it in the pipeline if needed.

The class currently provides a template, ``temp.pbs``, file to be
modified and used when submitting a job.

Code Examples
^^^^^^^^^^^^^

Submit 1 job
''''''''''''

.. code:: python

Submit multiple jobs
''''''''''''''''''''

.. code:: python

Get Job Info
''''''''''''

.. code:: python

:exclamation: Software Dependencies
-----------------------------------

Ensure that you have at least pbs version ``14.1.0``

Thanks
------

Thanks to [@jfeala](https://github.com/jfeala) for his work on Luigi's
SGEJobTask.
