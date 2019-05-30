Manager Documentation
=====================

The classes and functions in this module have been designed to help
manage existing and new projects using the Cookies module as well as the
different utilities found in the Tools module.

Why a manager?
--------------

This module is intended to mesh with a Flask user interface. \* Whenever
a new website is made the RepoManagement and WebManagement classes are
used. \* Whenever a new user is created in the Flask webpage, the
UserManagement class is used. \* Whenever an existing user creates a new
project, the ProjectManagement class is used.

However, this module does not have to be used to create a Flask webpage.
The full repository can be used for higher level organization, or
standalone projects can be made using the ProjectManagements
``_basic_project_`` flag.

The ``DataManagement`` class helps to tie everything together into a
pipeline.

Examples
--------

**Beware that this is under heavy development.**

Utilizing DataManagement to run a pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    import os
    from OrthoEvol.Manager import DataManagement

    DataManagement(pipeline="Ortho_CDS_1", start=True, new=True)

Utilizing DatabaseManagement to download databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

Notes
-----

Please view our detailed `BioSQL <biosqlreadme.html>`__ documentation and view
some of the static/config related
`files <https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager/config/>`__.
