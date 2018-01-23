Orthologs Documentation
=======================

This top level module includes submodules such as
`Align <https://github.com/datasnakes/OrthoEvolution/blob/master/OrthoEvol/Orthologs/Align/README.md>`__
(for aligning multi fasta files),
`Phylogenetics <https://github.com/datasnakes/OrthoEvolution/blob/master/OrthoEvol/Orthologs/Phylogenetics/README.md>`__
(for analyzing multiple sequence alignments), `BioSQL <>`__ (for
database creation),
`Blast <https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Orthologs/Blast>`__
(includes tools for using NCBI's blastn command line), and
`Genbank <https://github.com/datasnakes/OrthoEvolution/blob/master/OrthoEvol/Orthologs/Genbank/README.md>`__.
(for tools to extract features from genbank files).

Usage & Examples
----------------

These classes are optimized to be used together (very little work to do
that), but can also be used as standalone classes/methods.

This is a simple example of using all of the ``Orthologs`` submodules
together.

.. code:: python

    from OrthoEvol.Orthologs.Blast import OrthoBlastN
    from OrthoEvol.Orthologs.Align import ClustalO
    from OrthoEvol.Orthologs.Phlogenetics import ETE3PAML

Software Dependencies
---------------------

Ensure that the following software is installed and in your path:
Clustal omega, NCBI Blast+ 2.6.0 or greater, PAML, PhyML, Phylip,
IQTREE, Mafft, Prank, Clustalw, Guidance2 & Pal2Nal

If you are a sudo user, you may use the script we've provided,
`install.sh <https://github.com/datasnakes/OrthoEvolution/blob/master/OrthoEvol/Orthologs/install.sh>`__.

Using install.sh on Debian/Ubuntu:
----------------------------------

.. code:: bash

    # Change to the directory of the file.
    cd
    chmod +x install.sh
    ./sudo-install.sh

View Readmes for Orthologs Modules
------------------------------------
- `OrthoEvol.Orthologs.Align` - `Align <alignreadme.html>`__ documentation.
- `OrthoEvol.Orthologs.GenBank` - `GenBank <genbankreadme.html>`__ documentation.
- `OrthoEvol.Orthologs.Blast` - `Blast <blastreadme.html>`__ documentation.
- `OrthoEvol.Orthologs.Phylogenetics` - `Phylogenetics <phylogeneticsreadme.html>`__ documentation.
