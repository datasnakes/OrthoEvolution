.. image:: https://github.com/datasnakes/OrthoEvolution/actions/workflows/ci.yml/badge.svg
    :target: https://github.com/datasnakes/OrthoEvolution/actions/workflows/ci.yml

.. image:: https://badge.fury.io/py/OrthoEvol.svg
   :target: https://badge.fury.io/py/OrthoEvol

.. image:: https://readthedocs.org/projects/orthoevolution/badge/?version=latest
   :target: http://orthoevolution.readthedocs.io/en/latest/?badge=latest

.. image:: https://codecov.io/gh/datasnakes/OrthoEvolution/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/datasnakes/OrthoEvolution

.. image:: https://badgen.net/github/last-commit/datasnakes/OrthoEvolution
  :target: https://github.com/datasnakes/OrthoEvolution/commits/main




OrthoEvolution
====================
OrthoEvolution is an **easy to use** and comprehensive Python package which aids in the **analysis and
visualization of comparative evolutionary genetics** related projects such as the **inference of orthologs**.

**Current Version:** 1.0.0b2

Overview
--------------------------
This package focuses on **inferring orthologs** using NCBI's blast, various sequence alignment strategies,
and phylogenetics analyses including PAML, PhyML, ete3, and more tools.

Ultimately, the goal of this project is to create a **reusable pipeline** for the
inference of orthologs in order to ensure reproducibility of data as well as improve
the management and analysis of (what can be) large datasets.  The Cookies, Manager, Pipeline,
and Tools modules act as a framework for our workflow, while the Orthologs
module provides access to specific functions for our various ortholog inference projects.

View our `read the docs <http://orthoevolution.readthedocs.io/en/latest/>`__ and feel free to also
read `this related paper <https://www.frontiersin.org/articles/10.3389/fnhum.2014.00283/full>`__ to gain
more insight into this project/python package.


Installation
----------------
View the below methods for installing this package. Python 3.9 or higher is required.

PyPI
~~~~~~~~~~~~~~~~
.. code:: bash

    pip install --upgrade pip
    pip install OrthoEvol

GitHub
~~~~~~~~~~~~~~~~
.. code:: bash

    git clone https://github.com/datasnakes/OrthoEvolution.git
    cd OrthoEvolution
    pip install --upgrade pip
    pip install .

Development Code
~~~~~~~~~~~~~~~~
**WARNING** : This code is actively under development and may not be reliable. Please create an `issue <https://github.com/datasnakes/OrthoEvolution/issues>`_ for questions about development.

.. code:: bash

    git clone -b dev https://github.com/datasnakes/OrthoEvolution.git
    cd OrthoEvolution
    pip install --upgrade pip
    pip install .

Examples
----------------
Please view the `examples directory <https://github.com/datasnakes/OrthoEvolution/tree/main/examples>`__ for working examples and scripts demonstrating how to utilize this package.

The examples include:
- Standalone scripts for common workflows
- Example data files
- GUI implementations (Tkinter and PyWebView)
- Pipeline demonstrations


Running a pre-configured local blast
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: python

    from OrthoEvol.Orthologs.Blast import OrthoBlastN

    # Use an existing list of gpcr genes
    gpcr_blastn = OrthoBlastN(project="orthology-gpcr", method=1,
                             save_data=True, acc_file="gpcr.csv", 
                             copy_from_package=True)

    # Run blast
    gpcr_blastn.run()


Simple project creation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: python

    from OrthoEvol.Manager.management import ProjectManagement

    ProjectManagement(repo="test-repo", user=None,
                      project="test-project",
                      research=None,
                      research_type='comparative_genetics',
                      new_repo=False, new_user=False, new_project=True,
                      new_research=False)

Simple blast database downloading
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code:: python

    from OrthoEvol.Tools.ftp import NcbiFTPClient

    ncbiftp = NcbiFTPClient(email='somebody@gmail.com')
    ncbiftp.getblastdb(database_name='refseq_rna', v5=True)

Creating projects and databases dynamically 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        'cmd': f'python3 {os.path.join(os.getcwd(), job_name + ".py")}',
        'email': 'n/a'
    }

    # Save the updated configuration to a YAML file
    config_file_path = project_manager.user_log / Path("upload_config.yml")
    with open(str(config_file_path), 'w') as config_file:
        yaml.dump(db_config, config_file, default_flow_style=False)

    # Initialize database dispatcher and execute dispatch functions
    db_dispatcher = DatabaseDispatcher(config_file_path, project_manager)
    db_dispatcher.dispatch(db_dispatcher.strategies, db_dispatcher.dispatcher, db_dispatcher.configuration)


Tests
----------------
To run tests, first install the test dependencies:

.. code:: bash

    pip install pytest pytest-cov

Then run the test suite:

.. code:: bash

    pytest tests

Contributors
----------------
This package was created by the Datasnakes.

-  Rob Gilmore \| Github: `@grabear <https://github.com/grabear>`__ \|
   `✉ <mailto:robgilmore127@gmail.com>`__
-  Shaurita Hutchins \| Github:
   `@sdhutchins <https://github.com/sdhutchins>`__ \|
   `✉ <mailto:sdhutchins@outlook.com>`__

If you would like to contribute to this package, install the package in development mode:

.. code:: bash

    pip install -e .

Check out our `contributing guidelines <https://github.com/datasnakes/OrthoEvolution/blob/main/CONTRIBUTING.rst>`__ for more information.


License
----------------
This project is licensed under the MIT License. See the `LICENSE <https://github.com/datasnakes/OrthoEvolution/blob/main/LICENSE>`__ file for details.

Citations
----------------
We're thankful to have a resource such as
`Biopython <http://biopython.org/wiki/Biopython>`__, which inspired this
package.

*Cock, P.J.A. et al. Biopython: freely available Python tools for
computational molecular biology and bioinformatics. Bioinformatics 2009
Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163
pmid:19304878*
