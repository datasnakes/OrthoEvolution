Blast Documentation
===================

This module uses `NCBI's standalone
blast <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__
to generate blastn results. The results are parsed for the best hit,
which are used to get accession numbers.

What is BLAST?
--------------

Per NCBI, the `Basic Local Alignment Search Tool
(BLAST) <https://blast.ncbi.nlm.nih.gov/Blast.cgi>`__ finds regions of
local similarity between sequences. The program compares nucleotide or
protein sequences to sequence databases and calculates the statistical
significance of matches. BLAST can be used to infer functional and
evolutionary relationships between sequences as well as help identify
members of gene families.

We use NCBI's blastn task to generate a best hit in order to infer
orthology which is under the umbrella of comparative genetics/genomics
Comparative genetics/genomics is a field of biological research in which
the genome sequences of different species — human, mouse, and a wide
variety of other organisms from bacteria to chimpanzees — are compared.

Using this package, we compared these
`genes <http://www.guidetopharmacology.org/targets.jsp>`__ of interest
across a group of
`species <ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/>`__.

How do we configure and run blast?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running blast is the most complex aspect of this package, but we've
found a way to simplify the **automation of blasting** while also
**limiting blast searches by organism**.

Before you use this function, you need ``NCBI Blast+`` must be installed
and in your path. Download the latest standalone blast executables from
`here <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`__.

The story of 2 blast methods - Seqids vs Windowmasking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Seqids
^^^^^^

Windowmasking
^^^^^^^^^^^^^

We have perfected the method of using a windowmasker file for each
taxonomy id of the organisms that we are analyzing. The blastn
executable can filter a query sequence using the windowmasker data
files. This option can be used to mask interspersed repeats that may
lead to spurious matches. The windowmasker data files should be
downloaded from the NCBI FTP site.

For information on how to set up a window masker database, read our
`setup tutorial <window_masker_setup.md>`__.

On a command line, the windowmasker function would look as such:

.. code:: bash

    blastn -query input -db database -window_masker_taxid 9606 -out results.txt

That requires you to have a ``WINDOW_MASKER_PATH`` variable in your
environment variables.

In python:

.. code:: python

In addition to using windowmasker data files, we also use a specifically
formatted ``accession file`` with our headers as ``Tier``, ``Gene``,
``Organism`` to store blast output and input. This allows for
distinguishing genes by families or features. The ``Tier`` header can be
omitted, but the other headers are requirements.

The Accession numbers are stored in a .csv file. The following table is
an example of how we format our blast input file.

+------------+----------+-----------------+-------------------+-----------------+----------------------+
| Tier       | Gene     | Homo\_sapiens   | Macaca\_mulatta   | Mus\_musculus   | Rattus\_norvegicus   |
+============+==========+=================+===================+=================+======================+
| 1          | ADRA1A   | NM\_000680.3    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| 2          | ADRA1B   | NM\_000679.3    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| 3          | ADRA1D   | NM\_000678.3    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| 4          | ADRA2A   | NM\_000681.3    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| Good       | ADRA2B   | NM\_000682.6    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| Bad        | CHRM1    | NM\_000738.2    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| Ugly       | CHRM2    | NM\_000739.2    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| Other      | CHRM3    | NM\_000740.2    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| GPCR       | CHRM5    | NM\_012125.3    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+
| Isoforms   | CNR1     | NM\_016083.4    |                   |                 |                      |
+------------+----------+-----------------+-------------------+-----------------+----------------------+

The .csv file requires some manual configuration, and, while tedious, it
is also currently fundamental for the API.

Below we have defined the headers:

-  **Tier**: The target genes need a ranking or categorization based on
   the experiment. These can be user defined or a preset tier system can
   be used. In the future the different tiers will allow the user to
   control the order that each gene is processed.

-  **Gene**: The genes are HGNC aliases for the target genes of
   interest. In the future we will be able to process the HGNC .csv file
   to further automate the creation of this template file.

-  **Query**: The query organism is placed into the 3rd column of the
   .csv file. In the example Homo sapiens is used. Each taxa is a string
   in the format of "*Genus\_species*". The query organism also has to
   have accession numbers for each gene. It is therefore highly
   important to pick a well annotated species for accurate analysis.

Examples
--------

The main class to use is ``OrthoBlastN`` in order to run blast. In order
to run ``OrthoBlastN`` without using our database management features,
``BLASTDB`` and ``WINDOW_MASKER_PATH`` paths must be set.

Performing Blast & Post-Blast Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Orthologs.Blast import OrthoBlastN
    import os

    # Create a blast configuration dictionary
    blast_cfg = {
                  "taxon_file": None,
                  "go_list": None,
                  "post_blast": True,
                  "template": None,
                  "save_data": True,
                  "copy_from_package": True,
                  "MAF": 'MAFV3.2.csv'
                   }


    path = os.getcwd()
    myblast = OrthoBlastN(proj_mana=None, project="blast-test", project_path=path, **blast_config)

    # If you want to immediately start blasting, set auto_start to True
    myblast.blast_config(myblast.blast_human, 'Homo_sapiens', auto_start=False)

Making the API available with Accession data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*TODO: This is unfinished.*

.. code:: python

    from OrthoEvol.Orthologs.CompGenetics import CompGenAnalysis

