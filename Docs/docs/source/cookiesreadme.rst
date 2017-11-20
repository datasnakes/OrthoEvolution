Cookies Documentation
=====================

For this project/package, we recommend using cookiecutter (along with
Flask) to set up your directory if you intend to create a web
app/interface for your project.

Cookies makes it very easy to do this.

Learn more about the
`cookiecutter <https://github.com/audreyr/cookiecutter>`__ package.

Examples
--------

The `Manager
module <https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager>`__
uses the *CookBook* and *Oven* classes as a primary means of
functioning.

Simple Implementation
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Cookies import Oven

    Kitchen = Oven(repo="repo", user="user", project="project", output_dir="project_path")
    Pantry = Kitchen.Ingredients
    Kitchen.bake_the_*()
