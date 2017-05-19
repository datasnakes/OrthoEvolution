# Used:
# https://github.com/pypa/sampleproject/blob/master/setup.py
# https://github.com/biopython/biopython/blob/master/setup.py
# TODO-ROB: 1.  Generate a Template BioSQL database using sqlite3

from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import sys

sys.stderr = open('err.txt', 'w')

home = path.abspath(path.dirname(__file__))
name = 'Datasnakes-Scripts'

PACKAGES = [
    'Align',
    'Orthologs',
    'CompGenetics',
    'Blast',
    'BioSQL',
    'Manager',
    'Cookies',
    'Docs',
    'Phylogenetics'
    #'Tools',
]

# Get the long description from the README file
def readme():
    with open(path.join(home, 'README.rst'), encoding='utf-8') as f:
        return f.read()

# Set up an initial configuration file
# //TODO-ROB: Set up a reconfiguration script for manual movements;  Make the script append to 'init.yaml'

# Setup the package
setup(
    name=name,
    author = 'Datasnakes',
    description="A project that will help to analyze orthologous genes.",
    version='0.1.0',
    long_description=readme(),
    url='https://github.com/datasnakes/Datasnakes-Scripts',
    license='?',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Programming Language :: Python :: 3',
        'Operating System :: Unix',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5'
    ],
    packages=find_packages(exclude=['Docs', 'Tools', 'Archive', 'Examples']),
    install_requires=[],
    include_package_data=True,
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose']
)

