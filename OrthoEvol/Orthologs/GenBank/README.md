# Genbank Documentation

Retrieve genbank files and extract specific features sucha as `cds` or `aa`. Also,
write the features to text files.

If you haven't worked with genbank files before or are unfamiliar with what they
look like, view a [sample genbank record](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html).

## Usage

The main class is `GenBank`.

### Notable File Formats

|Extension|	Meaning   |Notes                                           |
|---------|-----------|------------------------------------------------|
|fasta    | generic fasta | Any generic fasta file. Other extensions can be fas, fa, seq, fsa |
|fna      | fasta nucleic acid | Used generically to specify nucleic acids. |
|ffn      | FASTA nucleotide of gene regions | Contains coding regions for a genome. |
|faa      | fasta amino acid | Contains amino acids. A multiple protein fasta file can have the more specific extension mpfa. |
|frn      | FASTA non-coding RNA | Contains non-coding RNA regions for a genome, in DNA alphabet e.g. tRNA, rRNA |


## Examples

### Perform Genbank Feature Extraction

``` python

```
