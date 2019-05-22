PAML Documentation
==================

PAML (Phylogenetic Analysis by Maximum Likelihood) is a package of
programs for phylogenetic analyses of DNA or protein sequences using
maximum likelihood and is maintained by Ziheng Yang.

Why ETE?
--------

ETE is python package for building, comparing, annotating, manipulating
and visualising trees. It provides a comprehensive API and a collection
of command line tools, including utilities to work with the NCBI
taxonomy tree.

Model Selection and Default Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It's important to note the default parameters for ``ETE3PAML`` are as
follows: ``model='M1'``, \`workdir=\`\`\`.

Usage & Examples
----------------

A simple implementation of ETE3PAML
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML

    paml = ETE3PAML(alignmentfile='.ffn', speciestree='.nw', workdir='')

    paml.run(pamlsrc='path/to/codeml/binary', output_folder=None)

Pruning a tree for use with ETE3PAML
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML

    paml = ETE3PAML(alignmentfile='HTR1A.ffn', speciestree='speciestree.nw', workdir='')

    # Input a list of orgnanisms or an organisms csv file with header as 'Organisms'
    paml.prune_tree(organisms='organisms.csv')

    paml.run(pamlsrc='path/to/codeml/binary', output_folder=None)

