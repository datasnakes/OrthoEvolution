CompGenetics (Comparative Genetics) Documentation
-------------------------
Perform comparative genetics bioinformatics studies on [genes](http://www.guidetopharmacology.org/targets.jsp)
of interest across a group of [species](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/multiprocessing/).


Usage
-----

The main classes under CompGenetics are `BLASTAnalysis` and `CompGenAnalysis`.

#### Code Examples

##### Performing Blast Analysis

``` python
from Orthologs.CompGenetics import BLASTAnalysis

# Take a look at the required and default parameters
# The default arguments are template, taxon_file, post_blast, save_data
BLASTAnalysis(self, repo, user, project, research, research_type,
              template=None, taxon_file=None, post_blast=False, save_data=True)


```
##### Performing Comparative Genetics Analysis

``` python
from Orthologs.CompGenetics import CompGenAnalysis

```

Tests
-----

Describe and show how to run the tests with code examples.

:exclamation: Notes
-------------------

Explain or list any notable information about the contents of this folder.

Other
-----
