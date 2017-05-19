## Datasnakes-Orthologs

This package is a collection of the scripts related to an Orthologs Project.

## Description

The subdirectories here are used as modules in the Python3 package that we are developing.

## Installation

Soon we'll be able to `pip install datasnakes-orthologs` via a command line.

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
from Orthologs import Align, BioSQL, Blast, CompGenetics, Phylogenetics, Genbank

```

## Contributors
* Rob Gilmore | Github: [@grabear](https://github.com/grabear) | [:email:](mailto:robgilmore127@gmail.com)
* Shaurita Hutchins | Github: [@sdhutchins](https://github.com/sdhutchins) | Twitter: [@MavenNBA](https://twitter.com/MavenNBA/) | [:email:](mailto:sdhutchins@outlook.com)


### Citation
We're so thankful to have a resource such as [Biopython](). They inspired this package.

*Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878*