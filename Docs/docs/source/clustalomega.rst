Align Documentation
----------------------------

This Align module aids in aligning multiple sequence fasta files, and in particular,
it has been designed to optimize aligning orthologous mammal sequences.  AQUA and
other alignment applications such as muscle or tcoffee may be added later. We've
found that `clustal omega <http://www.ebi.ac.uk/Tools/msa/clustalo/help/faq.html>`__
is best for the sample size we presently use which includes about 66 sequences
per mutli fasta file.


Usage
~~~~~~~~~~~~~~~

Clustal Omega is mainly used to align our the nucleotides sequences. AQUA and
other alignment applications such as muscle or tcoffee may be added later.
It can also be used fairly easily with multiprocessing.

Code Example
^^^^^^^^^^^^^^^^^^^^^

This is a quick example to use this class.

.. code:: python

    from Datasnakes.Orthologs import Align

    gene_list = ['HTR1A', 'CCR5', 'DRD4']

    for gene in gene_list:
        Align.ClustalO(gene + "_multifasta.ffn", gene + "_aligned.fasta", gene + ".log")



Format Examples
^^^^^^^^^^^^^^^^^^^^^

**Directory name:** ``cds``

**CDS input filename:** ``HTR1A_Homo_sapiens_cds_nucl.fasta``

**Fasta file format:**

``>H_sapiens``
``ATGGATGTGGTCAACAGCCTTCTTGTGAATGGAAGCAACATCACCCCTCCTTGCGAACTAGGCATCGAAAATGAGACGCTTTTCTGCTTGGATCAGCCCCATTCATCCAAAGAGTGGCAGCCAGCTGTGCAGATTCTCTTGTATTCCTTGATATTCCTGCTCAGCGTGCTGGGGAACACGCTGGTCATCACGGTGCTGATTCGAAACAAGAGGATGAGAACCGTCACCAACATCTTCCTGCTGTCCCTGGCCATCAGCGACCTCATGCTCTGCCTCTTCTGTATGCCGTTCAACCTCATCCCCAACCTGCTCAAGGATTTCATCTTCGGGAACGCTGTTTGCAAGACCACCACCTACTTCATGGGCACCTCTGTGAGTGTATCCACTTTTAATCTGGTCGCCATATCTCTGGAGCGATACGGTGCGATTTGCAAACCCTTACAGTCCCGGGTCTGGCAGACGAAATCCCACGCTTTGAAGGTGATTGCTGCTACCTGGTGCCTCTCCTTTACCATCATGACTCCGTACCCAATTTATAGCAACTTGGTGCCTTTTACCAAAAATAACAACCAGACCGCAAATATGTGCCGCTTTCTACTGCCAAATGATGTTATGCAGCAGTCCTGGCACACGTTCCTGTTACTCATCCTCTTTCTTATTCCTGGAATCGTAATGATGGTAGCATATGGATTAATCTCTTTGGAACTTTACCAAGGAATAAAATTTGATGCTAGCCAGAAGAAGTCTGCTAGAGAAAGGAAACCGAGCGCCGGCAGCAGCGGCCGATACGAGGACAGTGACGGGTGTTACCTGCAGAAGTCCAAGCACCCACGCAAGCTGGAGCTTCAGCGCCTGTCGCCCGGCGGCAGCGGCAGGGTCAACCGCATCAGGAGCAGCAGCTCCGCGGCCAACCTGATGGCCAAGAAGCGGGTGATCCGCATGCTCATGGTCATCGTGGTCCTCTTCTTCCTGTGCTGGATGCCCATCTTCAGCGTCAACGCCTGGCGGGCCTATGACACGGCCTCCGCCGAGCGCCACCTCTCGGGGACCCCCATTTCCTTCATCCTCCTGCTCTCCTACACCTCCTCCTGCGTCAACCCCATCATCTACTGCTTCATGAACAAACGGTTCCGCCTCGGCTTCATGGCCACCTTCCCCTGCTGCCCCAACCCTGGTCCCCCAGGGGCGAGAGGAGAGCCGGGAGAGGAGGAGGAAGGCAGGACCACCGGGGCCTCGCTGTCCAGGTACTCATACAGCCACATGAGCGCCTCTGCCCCGCCCCCGTGA``

**CDS output file name:** ``HTR1A_cds_nucl_aligned.fasta``


Analysis Parameters
^^^^^^^^^^^^^^^^^^^^^

``infile="MASTER_" + gene + "_CDS1.ffn",``

``outfile=gene + "_aligned_cds_nucl.fasta",``

``seqtype="DNA",``

``max_hmm_iterations=2,``

``infmt="fasta",``

``outfmt="fasta",``

``iterations=3,``

``verbose=True,``

``threads=8,``

``force=True,``

``log=gene + "_alignment.log"``

Notes
^^^^^^^^^^^^^^^^^^^^^
2 spaces should be after the species/common name in the fasta file.
