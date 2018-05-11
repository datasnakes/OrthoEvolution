Align Documentation
===================

This module aids in aligning multiple sequence fasta files, and in
particular, it has been designed to optimize aligning orthologous
mammalian sequences. We've found that `clustal
omega <http://www.ebi.ac.uk/Tools/msa/clustalo/help/faq.html>`__ is best
for the sample size we presently use which includes about 66 sequences
per mutli fasta file.

In the process of aligning our sequences, we also researched methods for
better curating those sequences. We've added
`Guidance2 <http://guidance.tau.ac.il/>`__ and
`Pal2Nal <http://www.bork.embl.de/pal2nal/>`__ command line wrappers to
help us to remove poor sequences (guidance2) and to prep sequences
better for PAML analysis (pal2nal).

Examples
--------

Clustal Omega is mainly used to align our the cds sequences. It's best
to use clustal omega with amino acid sequences.

Using the MultipleSequenceAlignment class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MultipleSequenceAlignment class was carefully designed to optimize
the alignment strategies for orthology inference and speed.

.. code:: python

    from OrthoEvol.Orthologs.Align import MultipleSequenceAlignment as MSA

    # User may specify a project and project path
    msa = MSA(project=None, project_path=os.getcwd())

    # From there, the user can select an alignment strategy (clustalo, guidance2, or pal2nal)

    fastafile = 'HTR1A.ffn'

    msa.clustalo(infile=fastafile, outfile='HTR1A_aligned_clustal.ffn', outfmt="fasta")

    msa.guidance2(seqFile=fastafile, msaProgram='MUSCLE', seqType='aa', dataset='MSA',
                  seqFilter=None, columnFilter=None, maskFilter=None)

    # AFTER aligning, the user may use pal2nal
    msa.pal2nal(aa_alignment='HTR1A_aligned_clustal.ffn', na_fasta=fastafile,
                output_type='paml', nogap=True, nomismatch=True, downstream='paml')

Running Clustal Omega
~~~~~~~~~~~~~~~~~~~~~

It's important to note that the default parameters for ``ClustalO`` are
as follows: ``seqtype="PROTEIN"``, ``infmt="fasta"``, ``outfmt="fasta"``

.. code:: python

    from OrthoEvol.Orthologs.Align import ClustalO

    gene_list = ['HTR1A', 'CCR5', 'DRD4']

    for gene in gene_list:
        ClustalO(infile=gene + "_multifasta.ffn", outfile=gene + "_aligned.fasta", logpath=gene + ".log")

Running the Pal2Nal command line wrapper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the out format is in clustal or aln format.

.. code:: python

    from OrthoEvol.Orthologs.Align import Pal2NalCommandline

    Pal2NalCommandline(pepaln='HTR1A_aligned.aln', nucfasta='HTR1A.ffn', output_file='HTR1A_aligned_pal2nal.aln')

Running the Guidance2 command line wrapper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The user can choose a multiple sequence alignment program (Options
include 'MAFFT', 'PRANK', 'CLUSTALW', 'MUSCLE') and a sequence type
(seqType) which can be 'aa', 'nt', or 'codon'.

.. code:: python

    from OrthoEvol.Orthologs.Align import Guidance2Commandline

    Guidance2Commandline(seqFile='HTR1A.ffn', msaProgram='MUSCLE', seqType='aa',
                         outDir='path/of/output/dir')
