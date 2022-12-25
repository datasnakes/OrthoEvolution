#!/usr/bin/env python
"""Invoke tasks."""

import os
import json
import shutil
import webbrowser

from invoke import task

HERE = os.path.abspath(os.path.dirname(__file__))

# Load the settings from cookiecutter.json
with open(os.path.join(HERE, 'cookiecutter.json'), 'r') as fp:
    COOKIECUTTER_SETTINGS = json.load(fp)

# Match default value of website_name from cookiecutter.json
COOKIE = os.path.join(HERE, COOKIECUTTER_SETTINGS['website_name'])

# Path to autoapp.py file
AUTOAPP = os.path.join(COOKIE, 'autoapp.py')

# Path to dev requirements file
REQUIREMENTS = os.path.join(COOKIE, 'requirements', 'dev.txt')

@task
def build(ctx):
    """Build the cookiecutter."""

    # Run cookiecutter with no input
    ctx.run(f'cookiecutter {HERE} --no-input')

@task
def clean(ctx):
    """Clean out generated cookiecutter."""

    # Remove the cookiecutter directory if it exists
    if os.path.exists(COOKIE):
        shutil.rmtree(COOKIE)
        print(f'Removed {COOKIE}')
    else:
        print('App directory does not exist. Skipping.')

def _run_flask_command(ctx, command):
    """Run a Flask command."""

    # Run the specified Flask command
    ctx.run(f'FLASK_APP={AUTOAPP} flask {command}', echo=True)

@task(pre=[clean, build])
def test(ctx):
    """Run lint commands and tests."""

    # Install dev requirements
    ctx.run(f'pip install -r {REQUIREMENTS} --ignore-installed', echo=True)

    # Change to the cookiecutter directory
    os.chdir(COOKIE)

    # Run lint and test commands
    _run_flask_command(ctx, 'lint')
    _run_flask_command(ctx, 'test')

@task
def readme(ctx, browse=False):
    """Convert the README to HTML."""

    # Convert the README to HTML
    ctx.run("rst2html.py README.rst > README.html")

    # Open the HTML file in a web browser if specified
    if browse:
        webbrowser.open_new_tab('README.html')
