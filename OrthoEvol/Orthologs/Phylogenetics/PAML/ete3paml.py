import os

import pandas as pd
from ete3 import EvolTree, Tree

from OrthoEvol.Tools.logit import LogIt


class ETE3PAML(object):
    """Integration of ETE3 for using PAML's codeml.

    M1 model is best for orthology inferences.
    """

    def __init__(self, infile, species_tree, workdir, pamlsrc=None):
        """Initialize main variables/files to be used.

        Ensure that you have the correct path to your codeml binary. It should
        be in the paml `/bin`.

        :param infile: The input fasta file.
        :type infile: [type]
        :param species_tree: [description]
        :type species_tree: [type]
        :param workdir: [description]
        :type workdir: [type]
        :param pamlsrc: [description], defaults to None
        :type pamlsrc: [type], optional
        """
        # Set up the logger
        self.paml_log = LogIt().default(logname="paml", logfile=None)

        self.infile = infile
        self.species_tree = species_tree
        self.workdir = workdir
        self.pamlsrc = pamlsrc

        if not self.pamlsrc:
            # If user does not specify a path, assume it is in path.
            self.pamlsrc = "codeml"

        # Import your species tree
        self._speciestree = Tree(self.species_tree, format=1)

        # Import alignment file as string
        self.aln_str = self._import_alignment()

    def _import_alignment(self):
        """Import alignment file as string."""
        with open(self.infile, 'r') as alignment_file:
            alignment_str = alignment_file.read()
        return alignment_str

    def prune_tree(self, organisms_list, organisms_file=None, column_header="Organisms"):
        """Prune branches for species not in the alignment file.

        Keep branches in the species tree for species in the alignment file
        Some species may not be present in the alignment file due to lack of
        matching with blast or simply the gene not being in the genome.
        """
        # If an organisms file is used, import and convert to list.
        if organisms_file:
            organisms_df = pd.read_csv(organisms_file)
            organisms_list = list(organisms_df[column_header])

        branches_to_keep = []
        try:
            for organism in organisms_list:
                if organism in self.aln_str:
                    branches_to_keep.append(organism)
                else:
                    self.paml_log.warning('No sequence for %s.' % organism)
    
                self._speciestree.prune(branches_to_keep, preserve_branch_length=True)
        except ValueError as e:
            self.paml_log.exception(e)
            
        else:
            # Write the tree to a file if not a ValueError
            temp_tree_path = os.path.join(self.workdir, 'temptree.nw')
            self._speciestree.write(outfile=temp_tree_path)

    def run(self, outfile, tree="temptree.nw", model="M1"):
        """Run PAML using ETE.

        The default model is M1 as it is best for orthology inference in
        our case. You can use models `M2`, `M0`, `M3`.
        """
        
        # Import the newick tree
        tree = EvolTree(tree)

        # Import the alignment
        tree.link_to_alignment(self.alignmentfile)

        tree.workdir = self.workdir

        # Set the binpath of the codeml binary
        tree.execpath = self.pamlsrc
        
        tree.run_model(model + '.' + outfile)  # Run the model M1 M2 M3 M0
