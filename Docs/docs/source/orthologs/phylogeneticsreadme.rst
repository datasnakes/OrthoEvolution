Phylogenetics Documentation
===========================

This documentation will provide information and guidelines about how we
use the Phylogenetics modules related to this package.

Overview
--------

Phylogenetics is best defined as the study of evolutionary relationships
among biological entities. In our case, those entities are species. We
are seeking to learn how mammals (more specifically primates) compare to
each other given a group of genes (GPCRs and addiction related).

PAML in particular is the most rigorous software for helping us to
understand the potentially significant differences in genes across
different mammalian species. From there, we can decide which genes we
will further study in cell culture projects or assays.

Dependencies
------------

It is critical to have these phylogenetic software installed and
avalailable on your path in order to use the Phylogenetics submodules or
you can take a look at our `external apps
repository <https://github.com/datasnakes/external-apps>`__ to help you
install these software on your machine.

Modules/Attributes available
----------------------------

.. code:: python

    from OrthoEvol.Orthologs import Phylogenetics

    # Find out what subclasses are available for use
    dir(Phylogenetics)

    Out[1]:
    ['AlignIO',
     'ETE3PAML',
     'FilteredTree',
     'IQTree',
     'IQTreeCommandline',
     'OrthologsWarning',
     'PhyML',
     'Phylip',
     'RelaxPhylip',
     'TreeViz',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__',
     'warnings']

Examples
--------

In the beginning stages of our project, we tested various phylogenetic
programs to see which worked well for us.

In this module, we include classes and ways to use PAML, Phylip, PhyML,
IQTREE, and Biopython's Bio.Phylo class.

Using PhyML with RelaxPhylip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    # Now you can import a class you want to utilize
    from OrthoEvol.Orthologs.Phylogenetics import PhyML, RelaxPhylip

    RelaxPhylip("HTR1A_aligned.fasta", "HTR1A_aligned.phy")

    # Generate a maximum likelihood tree from the phylip formatted alignment file.
    PhyML("HTR1A_aligned.phy")

View detailed `PhyML <phymlreadme.html>`__ documentation.


Using ETE3PAML
~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML

    paml = ETE3PAML(alignmentfile='.ffn', speciestree='.nw', workdir='')

    paml.run(pamlsrc='path/to/codeml/binary', output_folder=None)

View detailed `PAML <pamlreadme.html>`__ documentation.

Using the FilteredTree implementation of IQTree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics import FilteredTree

    FilteredTree(alignment, dataType='CODON', working_dir='path/of/working/directory')

View detailed `IQTree <iqtreereadme.html>`__ documentation.
