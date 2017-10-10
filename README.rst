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
An **easy to use** and comprehensive package which aids in the **analysis and visualization of comparative genomics** & related bioinformatics projects.
This package places emphasis on the **inference of orthologs** using NCBI's blast, aligning sequences,
and phylogenetics analyses including PAML, PhymL, ete3, and more tools.

The goal of this project was to create a simple, short, & effective pipeline (refseq accession to phyloanalysis) to infer orthologs.

Check out our `wiki docs <https://github.com/datasnakes/Datasnakes-Scripts/wiki>`__ and our `read the docs <http://datasnakes-scripts.readthedocs.io/en/master/>`__!


Installation
------------
``pip install Datasnakes-Scripts``

Alternatively, you can set install the package manually.

1. Download the zip file and unzip it or ``git clone https://github.com/datasnakes/Datasnakes-Scripts.git``
2. ``cd Datasnakes-Scripts``
3. ``pip install .``


Examples
---------

.. code:: python

    from Datasnakes.Orthologs import Align, BioSQL, Blast, CompGenetics, Phylogenetics, Genbank

Tests
------
To run tests, type ``nosetests Tests/`` in the Datasnakes-Orthologs directory.

Contributors
------------

-  Rob Gilmore \| Github: `@grabear <https://github.com/grabear>`__ \|
   `:email: <mailto:robgilmore127@gmail.com>`__
-  Shaurita Hutchins \| Github:
   `@sdhutchins <https://github.com/sdhutchins>`__ \| Twitter:
   `@MavenNBA <https://twitter.com/MavenNBA/>`__ \|
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
