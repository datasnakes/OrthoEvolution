Blast Documentation
=====================
This package uses [NCBI's standalone blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
to generate blastn results.  The results are parsed for the best hit,
which are used to get accession numbers.

The Accession numbers are stored in a .csv file.  The following table is a good example.  Take note of the headers.

Tier      |  Gene    |  Homo_sapiens  |  Macaca_mulatta  |  Mus_musculus  |  Rattus_norvegicus
----------|----------|----------------|------------------|----------------|-------------------
1         |  ADRA1A  |  NM_000680.3   |                  |                |
2         |  ADRA1B  |  NM_000679.3   |                  |                |
3         |  ADRA1D  |  NM_000678.3   |                  |                |
4         |  ADRA2A  |  NM_000681.3   |                  |                |
Good      |  ADRA2B  |  NM_000682.6   |                  |                |
Bad       |  CHRM1   |  NM_000738.2   |                  |                |
Ugly      |  CHRM2   |  NM_000739.2   |                  |                |
Other     |  CHRM3   |  NM_000740.2   |                  |                |
GPCR      |  CHRM5   |  NM_012125.3   |                  |                |
Isoforms  |  CNR1    |  NM_016083.4   |                  |                |

The .csv file requires some manual configuration, and, while tedious, it is
also currently fundamental for the API.

Below we have defined the headers:

* **Tier**:  The target genes need a ranking or categorization based on the
experiment.  These can be user defined or a preset tier system can be used.
In the future the different tiers will allow the user to control the order
that each gene is processed.
* **Gene**:  The genes are HGNC aliases for the target genes of interest.
In the future we will be able to process the HGNC .csv file to further
automate the creation of this template file.
* **Query**:  The query organism is placed into the 3rd column of the .csv
file.  In the example Homo sapiens is used.  Each taxa is a string in the
format of "_Genus\_species_".  The query organism also has to have
accession numbers for each gene.  It is therefore highly important to pick a
well annotated species for accurate analysis.

What is Comparative Genetics?
-----------------------------
Comparative genetics/genomics  is a field of biological research in which the
genome sequences of different species — human, mouse, and a wide variety of
other organisms from bacteria to chimpanzees — are compared

For our lab, we compare these [genes](http://www.guidetopharmacology.org/targets.jsp)
of interest across a group of [species](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/).

What is BLAST?
----------------
Per NCBI, the [Basic Local Alignment Search Tool (BLAST)](https://blast.ncbi.nlm.nih.gov/Blast.cgi) finds regions of local
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
# First use the Manager module to set up directories

from Datasnakes.Manager import ProjectManagement

# This is more pythonic with YAML loading
Management_config = {
  "new_repo": True,
  "new_user": True,
  "new_project": True,
  "new_database": True,
  "new_research": True,
  "new_app": False,
  "new_website": False,
  "database": 'Test-Database',
  "repo": 'Test-Repository',
  "user": 'grabear',
  "project": 'Test-Project',
  "research": 'Test-Research',
  "research_type": 'Test-Research-Type'
}
pm = ProjectManagement(**Management_config)

# Second use the Blast module to start a blast
# Note:  The current API is being used for development

from Datasnakes.Orthologs.Blast import CompGenBLASTn

# This is more pythonic with YAML loading
Blast_config = {
  "taxon_file": None,
  "go_list": None,
  "post_blast": True,
  "template": None,
  "save_data": True,
  "copy_from_package": True,
  "MAF": 'MAFV3.2.csv'
}
bl = CompGenBLASTn(proj_mana=pm, **Management_config, **Blast_config)
bl.blast_config(bl.blast_human, 'Homo_sapiens', auto_start=True)

```
##### Making the API available with Accession data
_TODO: This is unfinished._

``` python
from Datasnakes.Orthologs.CompGenetics import CompGenAnalysis

```
Tests
-----

Describe and show how to run the tests with code examples.

:exclamation: Notes
-------------------
- [ ] Explain Window Masker.
