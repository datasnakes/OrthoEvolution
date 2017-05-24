Align Documentation
-------------------------
This module aids in aligning multiple sequence fasta files, and in particular,
it has been designed to optimize aligning orthologous mammal sequences.  AQUA and
other alignment applications such as muscle or tcoffee may be added later. We've
found that [clustal omega](http://www.ebi.ac.uk/Tools/msa/clustalo/help/faq.html)
is best for the sample size we presently use which includes about 66 sequences
per mutli fasta file.

Usage
-----

Clustal Omega is mainly used to align our the nucleotides sequences. AQUA and
other alignment applications such as muscle or tcoffee may be added later.
It can also be used fairly easily with multiprocessing.

#### Code Example

This is a quick example to use this class.

``` python
from Orthologs.Align import ClustalO

gene_list = ['HTR1A', 'CCR5', 'DRD4']

for gene in gene_list:
    ClustalO(gene + "_multifasta.ffn", gene + "_aligned.fasta", gene + ".log")

```


:exclamation: Notes
-------------------

It's important to not that the default parameters are as follows:
`seqtype="DNA"`, `infmt="fasta"`, `outfmt="fasta"`

