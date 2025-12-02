PhyML Documentation
===================

PhyML is a phylogeny software based on the maximum-likelihood principle.
Early PhyML versions used a fast algorithm performing Nearest Neighbor
Interchanges (NNIs) to improve a reasonable starting tree topology.

Learn more about PhyML `here <http://www.atgc-montpellier.fr/>`__.

Default Parameters
------------------

The default datatype is ``'aa' (amino acid)``, but you may use ‘nt’ for
nucleotide.

Examples
--------

Running Phyml
~~~~~~~~~~~~~

.. code:: python

   from OrthoEvol.Orthologs.Phylogenetics.PhyML import PhyML

   htr1a = PhyML(infile='HTR1A.phy', datatype='aa')
   htr1a.run()

Running Phyml with our parallel module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   from OrthoEvol.Orthologs.Phylogenetics.PhyML import PhyML
   from OrthoEvol.Tools.parallel import Multiprocess

   files = ['HTR1A.phy', 'HTR1E.phy', 'MAOA.phy']

   def phyml(filename):
       phyml = PhyML(infile=filename, datatype='aa')
       phyml.run()

   if __name__ == '__main__':
       mp = Multiprocess()
       mp.map2function(phyml, files)

Notes
-----

This class is designed for PhyML `version
3.1 <http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip>`__.
