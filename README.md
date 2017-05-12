## Datasnakes-Scripts

This package is a collection of the scripts related to an Orthologs Project. 

## Description

The subdirectories here are used as modules in the Python3 package that we are developing.

## Installation

In development but working.  To test you'll want to create a virtual environment so that cleanup is easy.
Using _virtualenv_ with python3 insures that _python_ invokes py3.5 and _python3_ invokes py36.  Invoke python36.

```bash
$ mkdir dev
$ cd dev
$ virtualenv PackageTest --python=python3
$ source activate PackageTest
$ cd PackageTest
$ pip install cookiecutter
$ git clone -b RAG-Review http;//github.com/datasnakes/Datasnakes-Scripts
$ cd Datasnakes-Scripts
$ python3 tester.py
```
## Usage

After installation, you'll be able to easily import each module via:

```python
from Orthologs import biosql, blast, ftp, genbank, phylogenetics
```


## Contributors
* Rob Gilmore | Github: [@grabear](https://github.com/grabear) | Email: [:email:](mailto:robgilmore127@gmail.com)
* S. Hutchins | Github: [@sdhutchins](https://github.com/sdhutchins) | Twitter: [@MavenNBA](https://twitter.com/MavenNBA/) | Email: [:email:](mailto:sdhutchins@outlook.com)
