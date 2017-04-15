## Vallender-Labs-Scripts Directory

This repository is a fork of the scripts directory in the Orthologs Project.  This will help to isolate some of the files and scripts for easier management during development.
## Description

The subdirectories here are used as modules in the Python3 package that we are developing.
Here is a link to the project directory that this repository represents:  https://github.com/robear22890/Orthologs-Project/tree/Robs-Rework/lib/scripts

## Usage

After installation via pip we will be able to easily import each module via:

from lib.scripts import manager, biosql, blast, ftp, genbank, multiprocessing

## Tests

To test this first install the package.  It is currently in development so a pull request will be necessary.
$ git clone https://github.com/robear22890/Orthologs-Project/tree/Robs-Rework
$ cd Orthlogs-Project
$ pip install .

In a live Python3.6 Console:

'''Python
from lib.scripts import *
'''