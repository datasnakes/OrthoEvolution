from ete3 import EvolTree, Tree
from OrthoTools import formatlist
import pandas as pd


class ETE3PAML(object):
    """Integration of ETE3 for using PAML's codeml. For this project,
    we start my using the M1 model so it is the default.
    """

    def __init__(self, gene, paml_path,
                 workdir='data/paml-output/', model='M1'):
        """Improve docstrings here."""
        # First make sure that
        t = Tree('data/initial-data/species_tree.nw', format=1)
        orgsfile = pd.read_csv('data/initial-data/organisms.csv', header=None)

        # Create a list name/variable and use list()
        orgs = list(orgsfile[0])
        organismslist = formatlist(orgs)

        # Import alignment file as string
        alignment_file = open(
            'data/clustal-output/' +
            gene +
            '_Aligned/' +
            gene +
            '_aligned_cds_nucl.fasta',
            'r')
        alignment_str = alignment_file.read()
        alignment_file.close()

        # Keep the branches in the species tree for species in the alignment file
        # Some species may not be present in the alignment file
        branches2keep = []
        for organism in organismslist:
            if organism in alignment_str:
                # print('Yup.')
                branches2keep.append(organism)
            else:
                pass
                # print('Nope.') Make an error code in the log

        # Input a list of branches to keep on the base tree
        speciestree = t.prune(branches2keep, preserve_branch_length=True)

        # Run PAML
        # Import the newick tree
        tree = EvolTree(speciestree)

        # Import the alignment
        tree.link_to_alignment(
            'data/clustal-output/' +
            gene +
            '_Aligned/' +
            gene +
            '_aligned_cds_nucl.fasta')

        tree.workdir = workdir

        # Set the binpath of the codeml binary
        tree.execpath = paml_path

        # Run the codeml model
        tree.run_model(model + '.' + gene)
