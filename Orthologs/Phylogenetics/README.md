# Phylogenetics Documentation

This document will provide information and guidelines about how we use the Phylogenetics modules related to this package.

Usage
-----

In the beginning stages of our project, we tested various phylogenetic programs
to see which worked well for us.

In this module, we include classes and ways to use PAML, Phylip, PhyML, and
Biopython's Bio.Phylo class.


#### Code Example

This is a quick example to use this class.

``` python
import Orthologs
from Phylogenetics.PAML import ete3paml

```


:exclamation: Notes
-------------------

It's important to not that the default parameters are as follows:
 `workdir='data/paml-output/'`, `model='M1'`

