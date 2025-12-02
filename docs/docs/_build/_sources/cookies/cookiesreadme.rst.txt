Cookies Documentation
=====================

For this package, we recommend using cookiecutter (along with Flask or
Dash which is built with Flask) to set up your directory if you intend
to create a web app/interface for your project.

``Cookies`` makes it very easy to do this.

Learn more about the
`cookiecutter <https://github.com/audreyr/cookiecutter>`__ package.

Overview
--------

The Cookies module provides two main classes:

- **CookBook**: Manages cookiecutter template paths and configurations
- **Oven**: Deploys cookiecutter templates to create directory
  structures

The `Manager
module <https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager>`__
uses the *CookBook* and *Oven* classes as a primary means of
functioning.

Available Templates
-------------------

- Templates used when creating a full repository:

  - ``new_repository`` - Creates a full repository structure
  - ``new_user`` - Creates a user directory structure
  - ``new_project`` - Creates a project directory structure
  - ``new_research`` - Creates a research directory structure
  - ``new_database_repo`` - Creates database repository structure
  - ``new_app`` - Creates R-Shiny application structure
  - ``new_website`` - Creates Flask website structure

- Template for standalone projects:

  - ``new_basic_project`` - Creates a basic standalone project structure

Examples
--------

Simple Implementation
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from OrthoEvol.Cookies import Oven

   # Create an Oven instance
   Kitchen = Oven(repo="repo", user="username", project="project-name",
                  output_dir="path/to/project")

   # Access ingredients (parameters)
   Pantry = Kitchen.Ingredients

   # Create different directory structures
   Kitchen.bake_the_repo()      # Create repository
   Kitchen.bake_the_user()      # Create user directory
   Kitchen.bake_the_project()   # Create project directory

Full Repository Setup
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from OrthoEvol.Cookies import Oven
   from pathlib import Path
   import os

   # Set up paths and names
   home = os.getcwd()
   repo = "Development"
   user = "username"
   project = "MyProject"
   research = "GPCR"
   research_type = "comparative_genetics"

   # Create paths
   repo_path = Path(home) / Path(repo)
   user_path = repo_path / Path('users')
   project_path = user_path / Path(user) / Path('projects')
   research_path = project_path / Path(project)

   # Initialize Oven for full repository
   Full_Kitchen = Oven(repo=repo, user=user, project=project, 
                       basic_project=False, output_dir=home)

   # Create the full structure
   Full_Kitchen.bake_the_repo()
   Full_Kitchen.bake_the_user(cookie_jar=user_path)
   Full_Kitchen.bake_the_project(cookie_jar=project_path)
   Full_Kitchen.bake_the_research(research=research, 
                                  research_type=research_type, 
                                  cookie_jar=research_path)

Basic Standalone Project
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from OrthoEvol.Cookies import Oven

   # Initialize Oven for basic project
   Basic_Kitchen = Oven(project="MyProject", basic_project=True, 
                        output_dir=os.getcwd())

   # Create the basic project
   Basic_Kitchen.bake_the_project()

Available Methods
~~~~~~~~~~~~~~~~~

The ``Oven`` class provides the following methods:

- ``bake_the_repo(cookie_jar=None)`` - Create a repository structure
- ``bake_the_user(cookie_jar=None)`` - Create a user directory structure
- ``bake_the_project(cookie_jar=None)`` - Create a project directory
  structure
- ``bake_the_research(research_type, research, cookie_jar=None)`` -
  Create a research directory
- ``bake_the_db_repo(db_config_file, db_path, cookie_jar=None, archive_flag=False, delete=False)``
  - Create database repository
- ``bake_the_website(host, port, website_path, cookie_jar=None)`` -
  Create a Flask website
- ``bake_the_app(app, cookie_jar=None)`` - Create an R-Shiny app
  structure
