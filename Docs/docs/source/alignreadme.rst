Align Documentation
-------------------

This module aids in aligning multiple sequence fasta files, and in
particular, it has been designed to optimize aligning orthologous
mammalian sequences. We've found that `clustal
omega <http://www.ebi.ac.uk/Tools/msa/clustalo/help/faq.html>`__ is best
for the sample size we presently use which includes about 66 sequences
per mutli fasta file.

In the process of aligning our sequences, we also researched methods for
better curating those sequences. We've added `Guidance2 <>`__ and
`Pal2Nal <>`__ command line wrappers to help us to remove poor sequences
(guidance) and to prep sequences better for PAML analysis (pal2nal).

Usage
-----

Clustal Omega is mainly used to align our the cds sequences. It's best
to use clustal omega with amino acid sequences.

Code Example
^^^^^^^^^^^^

This is a quick example to use the `ClustalO <>`__ class.

.. code:: python

    from Datasnakes.Orthologs.Align import ClustalO

    gene_list = ['HTR1A', 'CCR5', 'DRD4']

    for gene in gene_list:
        ClustalO(gene + "_multifasta.ffn", gene + "_aligned.fasta", gene + ".log")

:exclamation: Notes
-------------------

It's important to not that the default parameters are as follows:
``seqtype="PROTEIN"``, ``infmt="fasta"``, ``outfmt="fasta"``
