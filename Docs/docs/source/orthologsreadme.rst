Orthologs Documentation
=======================

This top level module includes submodules such as
`Align <https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Align/README.md>`__
(for aligning multi fasta files),
`Phylogenetics <https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Phylogenetics/README.md>`__
(for analyzing multiple sequence alignments), `BioSQL <>`__ (for
database creation),
`Blast <https://github.com/datasnakes/Datasnakes-Scripts/tree/master/Datasnakes/Orthologs/Blast>`__
(includes tools for using NCBI's blastn command line), and
`Genbank <https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Genbank/README.md>`__.
(for tools to extract features from genbank files).

Usage
-----

These classes are optimized to be used together (very little work to do
that), but can also be used as standalone classes/methods.

Simple Example
^^^^^^^^^^^^^^

This is a simple example of using some of the modules.

.. code:: python

    from Datasnakes.Orthologs import Phylogenetics as Phylo

Software Dependencies
---------------------

Ensure that the following software is installed and in your path: 1.
Clustal omega 2. NCBI Standalone Blast 3. PAML 4. PhyML 5. Phylip 6.
IQTREE 7. Mafft 8. Prank 9. Clustalw 10. Guidance2 11. Pal2Nal

If you are a sudo user, you may use the script we've provided,
`install.sh <https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/install.sh>`__.

Using ``install.sh`` on Debian/Ubuntu:

.. code:: bash

    # Change to the directory of the file.
    cd
    chmod +x install.sh
    ./sudo-install.sh

