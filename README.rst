.. image:: https://travis-ci.org/datasnakes/Datasnakes-Scripts.svg?branch=SDH-Review
    :target: https://travis-ci.org/datasnakes/Datasnakes-Scripts

Datasnakes-Orthologs
--------------------

This package is a collection of the scripts related to an Orthologs
Project.

Description
-----------

The subdirectories here are used as modules in the Python3 package that
we are developing.

Dependencies
------------

Currently, this package works with python 3.4 and upward.

Installation
------------

Soon we'll be able to ``pip install datasnakes-orthologs`` via a command
line.

In development but working. To test you'll want to create a virtual
environment so that cleanup is easy. Using *virtualenv* with python3
insures that *python* invokes py3.5 and *python3* invokes py36. Invoke
python36.

.. code:: bash

    $ mkdir dev
    $ cd dev
    $ virtualenv PackageTest --python=python3
    $ source activate PackageTest
    $ cd PackageTest
    $ pip install cookiecutter
    $ git clone -b RAG-Review http://github.com/datasnakes/Datasnakes-Scripts
    $ cd Datasnakes-Scripts
    $ python3 tester.py

Usage
-----

After installation, you'll be able to easily import each module via:

.. code:: python

    from Orthologs import Align, BioSQL, Blast, CompGenetics, Phylogenetics, Genbank

Contributors
------------

-  Rob Gilmore \| Github: `@grabear <https://github.com/grabear>`__ \|
   `:email: <mailto:robgilmore127@gmail.com>`__
-  Shaurita Hutchins \| Github:
   `@sdhutchins <https://github.com/sdhutchins>`__ \| Twitter:
   `@MavenNBA <https://twitter.com/MavenNBA/>`__ \|
   `:email: <mailto:sdhutchins@outlook.com>`__

Citation
~~~~~~~~

We're so thankful to have a resource such as
`Biopython <http://biopython.org/wiki/Biopython>`__. They inspired this
package.

*Cock, P.J.A. et al. Biopython: freely available Python tools for
computational molecular biology and bioinformatics. Bioinformatics 2009
Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163
pmid:19304878*


.. image:: https://api.codacy.com/project/badge/Grade/9a4ce39423ed4458a0c7fa3610c81ba2
   :alt: Codacy Badge
   :target: https://www.codacy.com/app/sdhutchins/Datasnakes-Scripts?utm_source=github.com&utm_medium=referral&utm_content=datasnakes/Datasnakes-Scripts&utm_campaign=badger