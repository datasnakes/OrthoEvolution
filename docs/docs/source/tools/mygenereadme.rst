MyGene Documentation
====================

Our ``MyGene`` class is a wrapper around BioThings' MyGene.info.
`MyGene.info <http://mygene.info>`__ provides simple-to-use REST web
services to query/retrieve gene annotation data.

Currently, our ``MyGene`` class does not allow any additional fields or
species (than the default), but more features will be added in the near
future.

Examples
--------

Use Blast Master Accession File output with MyGene
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Manager.config import templates

    infile = pkg_resources.resource_filename(templates.__name__, 'test_blast.csv')

    MyGene(infile=infile, outfile='mygene_output.csv')

    # Query mygene using your input file
    MyGene.query_mygene()
