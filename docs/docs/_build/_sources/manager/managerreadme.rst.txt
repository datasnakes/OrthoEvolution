Manager Documentation
=====================

The classes and functions in this module have been designed to help
manage existing and new projects using the Cookies module as well as the
different utilities found in the Tools module.

Why a manager?
--------------

This module is intended to mesh with a Flask user interface. While the
Flask server/client interface is not currently set up, we are developing
the package so that it will be easier to implement.

- Whenever a new website is made the RepoManagement and WebManagement
  classes are used.

  - Whenever a new user is created in the Flask webpage, the
    UserManagement class is used.
  - Whenever an existing user creates a new project, the
    ProjectManagement class is used.

This module does not have to be used to create a Flask webpage. The full
repository can be used for higher level organization, or standalone
projects can be made using the ProjectManagements ``_basic_project_``
flag.

Example
-------

Using DatabaseDispatcher to set up the proper databases using a ``YAML``
config file. This example is based on ``db_mana_test.py``:

.. code:: python

   from OrthoEvol.Manager.management import ProjectManagement
   from OrthoEvol.Manager.database_dispatcher import DatabaseDispatcher
   from OrthoEvol.Manager.config import yml
   from pkg_resources import resource_filename
   from pathlib import Path
   import yaml
   import getpass
   from datetime import datetime as d
   import os

   # Define job name
   job_name = "jobname"

   # Function to load configuration from YAML file
   def load_config(file_name):
       file_path = resource_filename(yml.__name__, file_name)
       with open(file_path, 'r') as file:
           return yaml.load(file, Loader=yaml.FullLoader)

   # Load project management configuration
   pm_config = load_config("initialize_new.yml")
   project_manager = ProjectManagement(**pm_config["Management_config"])

   # Load and update database management configuration
   db_config = load_config("databases.yml")
   db_config.update(pm_config)

   # Configure NCBI RefSeq release settings
   ncbi_config = db_config['Database_config']['Full']['NCBI']['NCBI_refseq_release']
   ncbi_config['upload_number'] = 12
   ncbi_config['pbs_dict'] = {
       'author': getpass.getuser(),
       'description': 'This is a default pbs job.',
       'date': d.now().strftime('%a %b %d %I:%M:%S %p %Y'),
       'proj_name': 'OrthoEvol',
       'select': '1',
       'memgb': '6gb',
       'cput': '72:00:00',
       'wt': '2000:00:00',
       'job_name': job_name,
       'outfile': job_name + '.o',
       'errfile': job_name + '.e',
       'script': job_name,
       'log_name': job_name,
       'pbsworkdir': os.getcwd(),
       'cmd': f'python3.6 {os.path.join(os.getcwd(), job_name + ".py")}',
       'email': 'n/a'
   }

   # Save the updated configuration to a YAML file
   config_file_path = project_manager.user_log / Path("upload_config.yml")
   with open(str(config_file_path), 'w') as config_file:
       yaml.dump(db_config, config_file, default_flow_style=False)

   # Initialize database dispatcher and execute dispatch functions
   db_dispatcher = DatabaseDispatcher(config_file_path, project_manager)
   db_dispatcher.dispatch(db_dispatcher.strategies, db_dispatcher.dispatcher, db_dispatcher.configuration)

Notes
-----

Please view our `BioSQL
documentation <https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager/BioSQL/README.md>`__
and view some of the static/config related
`files <https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager/config/>`__.
