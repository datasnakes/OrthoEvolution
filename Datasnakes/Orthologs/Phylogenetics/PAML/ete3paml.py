import os
from ete3 import EvolTree, Tree

from Datasnakes.Tools.utils import csvtolist


class ETE3PAML(object):
    """Integration of ETE3 for using PAML's codeml.

    M1 model is best for orthology inferences.
    """

    def __init__(self, inputfile, speciestree, workdir=''):
        """Initialize main variables/files to be used."""
        self.inputfile = inputfile
        self.speciestree = speciestree
        self.workdir = workdir

        # Import your species tree
        self._speciestree = Tree(self.speciestree, format=1)
        # TODO import organisms list

        # Import alignment file as string
        alignment_file = open(self.alignmentfile, 'r')
        alignment_str = alignment_file.read()
        self.aln_str = alignment_str
        alignment_file.close()

    def prune_tree(self, organismslist, organisms_file=None):
        """Prune branches for species not in the alignment file.

        Keep branches in the species tree for species in the alignment file
        Some species may not be present in the alignment file due to lack of
        matching with blast or simply the gene not being in the genome.
        """

        if organisms_file:
            organismslist = csvtolist(organisms_file)


        branches2keep = []
        for organism in organismslist:
            if organism in self.aln_str:
                branches2keep.append(organism)
            else:
                print('No sequence for %s.' % organism)

            self._speciestree.prune(branches2keep, preserve_branch_length=True)

            # Write the tree to a file
            self._speciestree.write(outfile=os.path.join(self.workdir,
                                                         'temptree.nw'))

    def run(self, pamlsrc, outfile, model='M1'):
        """Run PAML using ETE."""
        # Import the newick tree
        tree = EvolTree('temptree.nw')

        # Import the alignment
        tree.link_to_alignment(self.alignmentfile)

        tree.workdir = self.workdir

        # Set the binpath of the codeml binary
        tree.execpath = pamlsrc

        tree.run_model(model + '.' + outfile)  # Run the model M1 M2 M3 M0
