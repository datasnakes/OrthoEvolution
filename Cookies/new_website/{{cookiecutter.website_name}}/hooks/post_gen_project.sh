#!/usr/bin/env bash

export {{ cookiecutter.website_name | upper }}_SECRET='Genetics is cool for secrets right?'
export FLASK_APP="{{ cookiecutter.website_path }}autoapp.py"
export FLASK_DEBUG="1"

cd {{ cookiecutter.website_path }}
pip install -r requirements.txt
pip install -r requirements/dev.txt
bower install
flask run -h {{ cookiecutter.website_host }} -p {{ cookiecutter.website_port }}

"""Quickstart
----------

First, set your app's secret key as an environment variable. For example,
add the following to ``.bashrc`` or ``.bash_profile``.

.. code-block:: bash

    export {{cookiecutter.website_name | upper}}_SECRET='something-really-secret'

Before running shell commands, set the ``FLASK_APP`` and ``FLASK_DEBUG``
environment variables ::

    export FLASK_APP=/path/to/autoapp.py
    export FLASK_DEBUG=1

Then run the following commands to bootstrap your environment ::

    git clone https://github.com/{{cookiecutter.github_username}}/{{cookiecutter.website_name}}
    cd {{cookiecutter.website_name}}
    pip install -r requirements/dev.txt
    bower install
    flask run

You will see a pretty welcome screen.

Once you have installed your DBMS, run the following to create your app's
database tables and perform the initial migration ::

    flask db init
    flask db migrate
    flask db upgrade
    flask run
"""

#TODO-ROB when actually deploying make sure to change the settings below
"""

Deployment
----------

In your production environment, make sure the ``FLASK_DEBUG`` environment
variable is unset or is set to ``0``, so that ``ProdConfig`` is used.

"""