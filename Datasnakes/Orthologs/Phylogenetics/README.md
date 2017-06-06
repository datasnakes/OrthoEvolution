Phylogenetics Documentation
=============================
This document will provide information and guidelines about how we use the
Phylogenetics modules related to this package.

Phylogenetics is best defined as the study of evolutionary relationships among
biological entities. In our case, those entities are species. We are seeking to
learn how mammals (more specifically primates) compare to each other given a group
of genes (GPCRs and addiction related).

PAML in particular is the best software for helping us to understand the potentially
significant differences in genes across different mammalian species. From there, we can
decide which genes we will further study in cell culture projects or assays.

Usage
-----

In the beginning stages of our project, we tested various phylogenetic programs
to see which worked well for us.

In this module, we include classes and ways to use PAML, Phylip, PhyML, and
Biopython's Bio.Phylo class.


#### Code Example

This is a quick example to use this class.

``` python
from Orthologs import Phylogenetics

# Find out what subclasses are available for use
dir(Phylogenetics)

Out[1]:
['ETE3PAML',
 'OrthologsWarning',
 'PAML',
 'PamlTest',
 'PhyML',
 'Phylip',
 'PhyloTree',
 'PhymlTest',
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

# Now you can import a class you want to utilize
from Orthologs.Phylogenetics import PhyML, RelaxPhylip

RelaxPhylip("HTR1A_aligned.fasta", "HTR1A_aligned.phy")

PhyML("HTR1A_aligned.phy")


```


Default Parameters
-------------------

It's important to note the default parameters for `ETE3PAML` are as follows:
 `workdir='data/paml-output/'`, `model='M1'`

