Comparative Genetics Documentation and Blast Documentation
=====================
Use [NCBI's standalone blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
to generate blastn results.

What is Comparative Genetics?
-----------------------------
_TODO: Add a better description._

Perform comparative genetics bioinformatics studies on [genes](http://www.guidetopharmacology.org/targets.jsp)
of interest across a group of [species](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/multiprocessing/).

What is BLAST?
----------------
Per NCBI, the Basic Local Alignment Search Tool (BLAST) finds regions of local
similarity between sequences. The program compares nucleotide or protein
sequences to sequence databases and calculates the statistical significance of
matches. BLAST can be used to infer functional and evolutionary relationships
between sequences as well as help identify members of gene families.

NCBI's BLASTN programs search nucleotide databases using a nucleotide query.

Usage
-----

The main classes under CompGenetics are `CompGenFiles` and `CompGenObjects`.

#### Code Examples

##### Performing Blast Analysis

``` python
from Datasnakes.Orthologs.Blast import CompGenBLASTn

# Take a look at the required and default parameters
# The default arguments are template, taxon_file, post_blast, save_data
```
##### Performing Comparative Genetics Analysis

``` python
from Datasnakes.Orthologs.CompGenetics import CompGenAnalysis

```
Tests
-----

Describe and show how to run the tests with code examples.

:exclamation: Notes
-------------------

Explain or list any notable information about the contents of this folder.
