Comparative Genetics Documentation
----------------------------------

Perform comparative genetics bioinformatics studies on
`genes <http://www.guidetopharmacology.org/targets.jsp>`__ of interest
across a group of
`species <ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/multiprocessing/>`__.

Usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main classes under CompGenetics are ``BLASTAnalysis`` and
``CompGenAnalysis``.

Code Examples
^^^^^^^^^^^^^

Performing Blast Analysis
'''''''''''''''''''''''''

.. code:: python

    from Datasnakes.Orthologs.CompGenetics import BLASTAnalysis

    # Take a look at the required and default parameters
    # The default arguments are template, taxon_file, post_blast, save_data
    BLASTAnalysis(self, repo, user, project, research, research_type,
                  template=None, taxon_file=None, post_blast=False, save_data=True)

Performing Comparative Genetics Analysis
''''''''''''''''''''''''''''''''''''''''

.. code:: python

    from Orthologs.CompGenetics import CompGenAnalysis

