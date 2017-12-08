IQTree Documentation
====================

IQTree is a fast and effective stochastic algorithm to infer
phylogenetic trees by maximum likelihood.

Read more about iqtree at their `github
repository <https://github.com/Cibiv/IQ-TREE>`__.

Examples
--------

Using the commandline wrapper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, iqtree wil auto-detect a datatype. You can also specify a
datatype. You can select a datatype from the following list: ``BIN``,
``DNA``, ``AA``, ``NT2AA``, ``CODON``, ``MORPH``.

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics import IQTreeCommandline

    IQTreeCommandline(alignment='path/to/alignment/file')

Using the FilteredTree class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The FilteredTree class simply returns the best tree (which is determined
by IQTree's algorithm).

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics import FilteredTree

    FilteredTree(alignment, dataType='CODON', working_dir='path/of/working/directory')

Generating a Consensus Tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python


