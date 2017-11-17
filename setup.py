""" This is the setup.py script for setting up the package and fulfilling any
necessary requirements.

References:
https://github.com/pypa/sampleproject/blob/master/setup.py
https://github.com/biopython/biopython/blob/master/setup.py
http://python-packaging.readthedocs.io/en/latest/index.html
"""
# Modules Used
from setuptools import setup, find_packages
from codecs import open  # To use a consistent encoding
from os import path
import os
import sys
import pkg_resources
from importlib import import_module

## Save the standard error of the setup file. This can be removed soon.
#sys.stderr = open('err.txt', 'w')

# Set the home path of the setup script/package
home = path.abspath(path.dirname(__file__))
name = 'OrthoEvol'


def readme():
    """Get the long description from the README file."""
    with open(path.join(home, 'README.rst'), encoding='utf-8') as f:
        return f.read()

# Setup the package by adding information to these parameters


setup(
    name=name,
    author='Rob Gilmore & Shaurita Hutchins',
    author_email='datasnakes@gmail.com',
    description="This package aids in the analysis of orthologous genes.",
    version='0.9.0a2',
    long_description=readme(),
    url='https://github.com/datasnakes/OrthoEvolution',
    license='MIT',
    keywords='bioinformatics science evolution orthology psychiatry genetics',
    download_url='',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Programming Language :: Python :: 3',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Framework :: Flask',
        'Framework :: R-Shiny',
        'Framework :: Cookiecutter'
    ],
    # Packages will be automatically found if not in this list.
    packages=find_packages(exclude=['Docs', 'Examples', 'Tests']),
    include_package_data=True,
    entry_points={
        'console_scripts': [
                'orthoevol=OrthoEvol.Orthologs.command_line:main'
        ]
    },
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose']
)

# Set up the permissions for the BioSQL Perl scripts
scripts = import_module("OrthoEvol.Manager.BioSQL.biosql_repo.scripts")
biosql_scripts = pkg_resources.resource_filename(scripts.__name__, "")
for file in os.listdir(biosql_scripts):
    if '.pl' in file:
        script_path = os.path.join(biosql_scripts, file)
        os.chmod(script_path, mode=755)
