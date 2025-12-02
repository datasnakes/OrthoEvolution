# PAML Documentation

PAML (Phylogenetic Analysis by Maximum Likelihood) is a package of
programs for phylogenetic analyses of DNA or protein
sequences using maximum likelihood and is maintained by Ziheng Yang.


## Why ETE?

ETE is a python package for building, comparing, annotating, manipulating and visualising
trees. It provides a comprehensive API and a collection of command line tools including 
utilities to work with the NCBI taxonomy tree.

### Model Selection and Default Parameters

It's important to note the default parameters for `ETE3PAML` are as follows:
`model='M1'`, `workdir=''`.

## Usage & Examples

### A simple implementation of ETE3PAML

```python
from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML

paml = ETE3PAML(alignmentfile='.ffn', speciestree='tree.nw', workdir='', 
                pamlsrc='path/to/codeml/binary')

paml.run(output_folder=None)
```

### Pruning a tree for use with ETE3PAML

```python
from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML

paml = ETE3PAML(infile='HTR1A.ffn', species_tree='speciestree.nw', workdir='')

# Input a list of orgnanisms or an organisms csv file with header as 'Organisms'
paml.prune_tree(organisms='organisms.csv')

paml.run(pamlsrc='path/to/codeml/binary', output_folder=None)

```