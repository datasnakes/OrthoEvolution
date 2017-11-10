Cookies Documentation
=====================

For this project/package, we recommend using cookiecutter (along with
Flask) to set up your directory if you intend to create a web
app/interface for your project.

Cookies makes it very easy to do this.

Learn more about
`cookicutter <https://github.com/audreyr/cookiecutter>`__.

Examples
--------

The `Manager
module <https://github.com/datasnakes/Datasnakes-Scripts/tree/master/Datasnakes/Manager>`__
uses the *CookieRecipes* and *Oven* classes as a primary means of
functioning.

Here is a basic implementation:

.. code:: python

    from Datasnakes.Cookies import Oven

    Kitchen = Oven(repo="repo", user="user", project="project", output_dir="project_path")
    Pantry = Kitchen.Ingredients
    Kitchen.bake_the_*()
