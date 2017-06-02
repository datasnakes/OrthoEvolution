.. Datasnakes-Orthologs documentation master file, created by
   sphinx-quickstart on Wed May 24 14:34:01 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Datasnakes-Orthologs's documentation!
================================================

.. module:: Datasnakes

.. image:: https://travis-ci.org/datasnakes/Datasnakes-Scripts.svg?branch=release
   :target: https://travis-ci.org/datasnakes/Datasnakes-Scripts
.. image:: https://api.codacy.com/project/badge/Grade/9a4ce39423ed4458a0c7fa3610c81ba2
   :target: https://www.codacy.com/app/sdhutchins/Datasnakes-Scripts?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=datasnakes/Datasnakes-Scripts&amp;utm_campaign=Badge_Grade
.. image:: https://badges.gitter.im/gitterHQ/gitter.png
   :target: https://gitter.im/datasnakes/Lobby
.. image:: https://badge.fury.io/py/Datasnakes-Scripts.svg
   :target: https://badge.fury.io/py/Datasnakes-Scripts
.. image:: https://readthedocs.org/projects/datasnakes-scripts/badge/?version=master
   :target: http://datasnakes-scripts.readthedocs.io/en/master/
.. image:: https://img.shields.io/badge/under-development-orange.svg
   :target: https://github.com/datasnakes/Datasnakes-Scripts

The Datasnakes-Orthologs package is a collection of the modules which aids in
the analysis and visualization of orthologs-related bioinformatics projects.

Check out our `wiki docs <https://github.com/datasnakes/Datasnakes-Scripts/wiki>`__!

Installation
------------

For easy installation, ``pip install Datasnakes-Scripts``

For development mode installation:
1. ``git clone http://github.com/datasnakes/Datasnakes-Scripts``
2. ``cd Datasnakes-Scripts``
3. ``pip install -e .``

Dependencies
~~~~~~~~~~~~~

Currently, this package is optimized for ``python v3.4`` and upward. It's also
dependent upon ``biopython v1.68``.

Manual Installation
~~~~~~~~~~~~~~~~~~~~~

Alternatively, you can set install the package manually.

1. Download the zip file and unzip it or ``git clone https://github.com/datasnakes/Datasnakes-Scripts.git``
2. ``cd Datasnakes-Scripts``
3. ``pip install .``

Project Setup
~~~~~~~~~~~~~~~~
This package is in development but still working. To test you'll want to create a virtual
environment so that cleanup is easy. Using *virtualenv* with python3
insures that ``python`` invokes py3.5 and ``python3`` invokes py36. Invoke
python36.

.. code:: bash

    $ mkdir dev
    $ cd dev
    $ virtualenv PackageTest --python=python3
    $ source activate PackageTest
    $ git clone http://github.com/datasnakes/Datasnakes-Scripts
    $ cd Datasnakes-Scripts
    $ python3 setup.py install or pip install .

Usage
-----

After installation, you'll be able to easily import each module via:

.. code:: python

    from Datasnakes.Orthologs import Align, BioSQL, Blast, CompGenetics, Phylogenetics, Genbank

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

Contents
--------
.. toctree::
    :maxdepth: 2

    compgenetics
    clustalomega
    phylogenetics
    namingconventions


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`