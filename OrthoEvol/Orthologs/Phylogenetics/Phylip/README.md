# Phylip Documentation

PHYLIP (the PHYLogeny Inference Package) is a package of programs for inferring 
phylogenies (evolutionary trees).  Methods that are available in the package 
include parsimony, distance matrix, and likelihood methods, including 
bootstrapping and consensus trees. Data types that can be handled include 
molecular sequences, gene frequencies, restriction sites and fragments, 
distance matrices, and discrete characters.

Learn more about Phylip [here](http://evolution.genetics.washington.edu/phylip.html).

## Examples

### Running Phylip

```python
from OrthoEvol.Orthologs.Phylogenetics.Phylip import Phylip

htr1a = Phylip(infile='HTR1A.phy')

# Generate a distance matrix
htr1a.dnadist(outfile="htr1a_dist.txt")
```

### Running Phylip with our parallel module

```python

```
