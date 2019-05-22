OrthoEvolution
====================

.. image:: https://travis-ci.org/datasnakes/OrthoEvolution.svg?branch=master
    :target: https://travis-ci.org/datasnakes/OrthoEvolution

.. image:: https://api.codacy.com/project/badge/Grade/25062944794a4d14b5cab274a885ac27
   :target: https://www.codacy.com/app/datasnakes/OrthoEvolution?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=datasnakes/OrthoEvolution&amp;utm_campaign=Badge_Grade

.. image:: https://img.shields.io/badge/chat-on%20gitter-753A88.svg
   :target: https://gitter.im/datasnakes/OrthoEvolution

.. image:: https://badge.fury.io/py/OrthoEvol.svg
   :target: https://badge.fury.io/py/OrthoEvol

.. image:: https://readthedocs.org/projects/orthoevolution/badge/?version=latest
   :target: http://orthoevolution.readthedocs.io/en/latest/?badge=latest


OrthoEvolution is an **easy to use** and comprehensive python package which aids in the **analysis and
visualization of comparative evolutionary genetics** related projects such as the **inference of orthologs**.

Overview
--------------------------
This package focuses on **inferring orthologs** using NCBI's blast, various sequence alignment strategies,
and phylogenetics analyses including PAML, PhyML, ete3, and more tools.

Ultimately, the goal of this project is to create a **reusable pipeline** for the
inference of orthologs in order to ensure reproducibility of data as well as improve
the management and analysis of (what can be) large datasets.  The Cookies, Manager, Pipeline,
and Tools modules act as a framework for our workflow, while the Orthologs
module provides access to specific functions for our various ortholog inference projects.

View our `read the docs <http://orthoevolution.readthedocs.io/en/master/>`__ and feel free to also
read `this related paper <https://www.frontiersin.org/articles/10.3389/fnhum.2014.00283/full>`__ to gain
more insight into this project/python package.


Installation
----------------
View the below methods for installing this package. **Python3.5 and higher is required.**

PyPi
~~~~~~~~~~~~~~~~
``pip install OrthoEvol``

GitHub
~~~~~~~~~~~~~~~~
1. Download the zip file and unzip it or ``git clone https://github.com/datasnakes/OrthoEvolution.git``
2. ``cd OrthoEvolution``
3. ``pip install .``

Development Code
~~~~~~~~~~~~~~~~
**WARNING** : This code is actively under development and may not be reliable.  Please create an `issue <https://github.com/datasnakes/OrthoEvolution/issues>`_ for questions about development.

1. Download the zip file and unzip it or ``git clone -b dev-master https://github.com/datasnakes/OrthoEvolution.git``
2. ``cd OrthoEvolution``
3. ``pip install .``

Examples
----------------
View `detailed documentation <#detailed-documentation>`__ below.

.. code:: python

    import OrthoEvol

Tests
----------------
To run tests, type ``nosetests Tests/`` in the OrthoEvolution directory.

Contributors
----------------
This package was created by the Datasnakes.

-  Rob Gilmore \| Github: `@grabear <https://github.com/grabear>`__ \|
   `✉ <mailto:robgilmore127@gmail.com>`__
-  Shaurita Hutchins \| Github:
   `@sdhutchins <https://github.com/sdhutchins>`__ \|
   `✉ <mailto:sdhutchins@outlook.com>`__

If you would like to contribute to this package, install the package in development mode,
and check out our `contributing guidelines <https://github.com/datasnakes/OrthoEvolution/blob/master/CONTRIBUTING.rst>`__.


Citations
----------------
We're so thankful to have a resource such as
`Biopython <http://biopython.org/wiki/Biopython>`__. They inspired this
package.

*Cock, P.J.A. et al. Biopython: freely available Python tools for
computational molecular biology and bioinformatics. Bioinformatics 2009
Jun 1; 25(11) 1422-3 http://dx.doi.org/10.1093/bioinformatics/btp163
pmid:19304878*

License
----------------
`MIT <https://github.com/datasnakes/OrthoEvolution/blob/master/LICENSE>`_


.. toctree::
    :hidden:
    :maxdepth: 2

    index


Detailed Documentation
=========================

.. include:: docscontents.rst

Indices & Tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`