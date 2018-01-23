PhyML Documentation
===================

PhyML is a phylogeny software based on the maximum-likelihood principle.
Early PhyML versions used a fast algorithm performing Nearest Neighbor
Interchanges (NNIs) to improve a reasonable starting tree topology.

Learn more about PhyML `here <http://www.atgc-montpellier.fr/>`__.

Default Parameters
------------------

The default dataype is ``'aa' (amino acid)``, but you may use 'nt' for
nuclueotide.

Examples
--------

Running Phyml
~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML

    PhyML(phyml_input='path/to/phylip/multisequencealignment', datatype='aa')

Running Phyml with our parallel module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Orthologs.Phylogenetics.PAML import ETE3PAML

    PhyML(phyml_input='path/to/phylip/multisequencealignment', datatype='aa')

Notes
-----

This class is designed PhyML version 3.1.
