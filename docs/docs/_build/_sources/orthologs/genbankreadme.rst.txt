Genbank Documentation
=====================

Retrieve genbank files and extract specific features sucha as ``cds`` or
``aa``. Also, write the features to text files.

If you haven’t worked with genbank files before or are unfamiliar with
what they look like, view a `sample genbank
record <https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`__.

Usage
-----

The main class is ``GenBank``.

Notable File Formats
~~~~~~~~~~~~~~~~~~~~

+-----------+------------+-------------------------------------------------+
| Extension | Meaning    | Notes                                           |
+===========+============+=================================================+
| fasta     | generic    | Any generic fasta file. Other extensions can be |
|           | fasta      | fas, fa, seq, fsa                               |
+-----------+------------+-------------------------------------------------+
| fna       | fasta      | Used generically to specify nucleic acids.      |
|           | nucleic    |                                                 |
|           | acid       |                                                 |
+-----------+------------+-------------------------------------------------+
| ffn       | FASTA      | Contains coding regions for a genome.           |
|           | nucleotide |                                                 |
|           | of gene    |                                                 |
|           | regions    |                                                 |
+-----------+------------+-------------------------------------------------+
| faa       | fasta      | Contains amino acids. A multiple protein fasta  |
|           | amino acid | file can have the more specific extension mpfa. |
+-----------+------------+-------------------------------------------------+
| frn       | FASTA      | Contains non-coding RNA regions for a genome,   |
|           | non-coding | in DNA alphabet e.g. tRNA, rRNA                 |
|           | RNA        |                                                 |
+-----------+------------+-------------------------------------------------+

Examples
--------

Perform Genbank Feature Extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python
