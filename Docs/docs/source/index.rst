.. image:: https://travis-ci.org/datasnakes/Datasnakes-Scripts.svg?branch=master
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

Datasnakes-Scripts
--------------------
An **easy to use** and comprehensive package which aids in the **analysis and visualization of comparative genetics** & related bioinformatics projects.
The current implementation of this package places an emphasis on the **inference of orthologs** using NCBI's blast, various sequence alignment strategies,
and phylogenetics analyses including PAML, PhymL, ete3, and more tools.

The goal of this project is to create a pipeline framework for current (inference of orthologs) and future (RNA-seek) projects.  The Cookies, Manager, Pipeline,
and Tools modules act as a framework for our workflow, while the Orthologs module provides access to specific functions for our various ortholog inference projects.

Check out our `wiki docs <https://github.com/datasnakes/Datasnakes-Scripts/wiki>`__ and our `read the docs <http://datasnakes-scripts.readthedocs.io/en/master/>`__!


Installation
------------

PyPi
====================
``pip install Datasnakes-Scripts``

GitHub
===========
1. Download the zip file and unzip it or ``git clone https://github.com/datasnakes/Datasnakes-Scripts.git``
2. ``cd Datasnakes-Scripts``
3. ``pip install .``

Development Code
==================
**WARNING** : This code is actively under development and may not be reliable.  Please create an `issue <https://github.com/datasnakes/Datasnakes-Scripts/issues>`_ for questions about development.

1. Download the zip file and unzip it or ``git clone -b dev-master https://github.com/datasnakes/Datasnakes-Scripts.git``
2. ``cd Datasnakes-Scripts``
3. ``pip install .``

Examples
---------

.. code:: python

    import Datasnakes

A short overview of the package can be found in the Datasnakes module.
For more advanced documentation see the README.md files in each of the sub-modules.  A good place to start would be the
Orthologs/Blast module for understanding the bioinformatics workflow.

Tests
------
To run tests, type ``nosetests Tests/`` in the Datasnakes-Scripts directory.

Contributors
------------

-  Rob Gilmore \| Github: `@grabear <https://github.com/grabear>`__ \|
   `:email: <mailto:robgilmore127@gmail.com>`__
-  Shaurita Hutchins \| Github:
   `@sdhutchins <https://github.com/sdhutchins>`__ \|
   `:email: <mailto:sdhutchins@outlook.com>`__

If you would like to contribute to this package, install the package in development mode,
and check out our `contributing guidelines <https://github.com/datasnakes/Datasnakes-Scripts/blob/master/CONTRIBUTING.rst>`__.


Citation
----------

We're so thankful to have a resource such as
`Biopython <http://biopython.org/wiki/Biopython>`__. They inspired this
package.

*Cock, P.J.A. et al. Biopython: freely available Python tools for
computational molecular biology and bioinformatics. Bioinformatics 2009
Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163
pmid:19304878*

License
---------
`MIT <https://github.com/datasnakes/Datasnakes-Scripts/blob/master/LICENSE>`_

Contents
--------
.. toctree::
    :maxdepth: 3
    datasnakesreadme
    cookiesreadme
    managerreadme
    orthologsreadme
    toolsreadme

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`