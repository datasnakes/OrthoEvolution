import os
import pandas as pd
from ete3 import EvolTree, Tree

from OrthoEvol.utilities import FullUtilities

# Set up csv to list function
csvtolist = FullUtilities().csvtolist


class ETE3PAML(object):
    """Integration of ETE3 for using PAML's codeml.

    M1 model is best for orthology inferences.
    """

    def __init__(self, infile, species_tree, workdir, pamlsrc=None):
        """Initialize main variables/files to be used.
        
        Ensure that you have the correct path to your codeml binary. It should
        be in the paml `/bin`.

        :param infile: [description]
        :type infile: [type]
        :param species_tree: [description]
        :type species_tree: [type]
        :param workdir: [description]
        :type workdir: [type]
        :param pamlsrc: [description], defaults to None
        :type pamlsrc: [type], optional
        """
        self.infile = infile
        self.species_tree = species_tree
        self.workdir = workdir
        self.pamlsrc = pamlsrc
        self.temp_tree = None

        if not self.pamlsrc:
            # If user does not specify a path, assume it is in path.
            self.pamlsrc = "codeml"

        # Import your species tree
        self._speciestree = Tree(self.species_tree, format=1)
        # TODO import organisms list

        # Import alignment file as string
        alignment_file = open(self.infile, 'r')
        alignment_str = alignment_file.read()
        self.aln_str = alignment_str
        alignment_file.close()

    def prune_tree(self, organisms_list, organisms_file=None, column_header="Organisms"):
        """Prune branches for species not in the alignment file.

        Keep branches in the species tree for species in the alignment file
        Some species may not be present in the alignment file due to lack of
        matching with blast or simply the gene not being in the genome.
        """

        if organisms_file:
            og_df = pd.read_csv(organisms_file)
            organismslist = list(og_df[column_header])

        branches2keep = []
        for organism in organismslist:
            if organism in self.aln_str:
                branches2keep.append(organism)
            else:
                print('No sequence for %s.' % organism)

            self._speciestree.prune(branches2keep, preserve_branch_length=True)

            # Write the tree to a file
            temp_tree_path = os.path.join(self.workdir, 'temptree.nw')
            self.temp_tree = 'temptree.nw'
            self._speciestree.write(outfile=temp_tree_path)

    def run(self, outfile, tree, model='M1'):
        """Run PAML using ETE.

        The default model is M1 as it is best for orthology inference in
        our case. You can use models `M2`, `M0`, `M3`.
        """
        
        # Import the newick tree
        tree = EvolTree(self.temp_tree)

        # Import the alignment
        tree.link_to_alignment(self.alignmentfile)

        tree.workdir = self.workdir

        # Set the binpath of the codeml binary
        tree.execpath = self.pamlsrc

        tree.run_model(model + '.' + outfile)  # Run the model M1 M2 M3 M0
