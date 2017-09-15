from ete3 import EvolTree, Tree
import pandas as pd
from Datasnakes.Manager.utils import FormatList


class ETE3PAML(object):
    """Integration of ETE3 for using PAML's codeml.

    M1 model is best for orthology inferences.
    """

    def __init__(self, gene, paml_path, workdir='data/paml-output/',
                 model='M1'):
        """Improve docstrings here."""
        # Import your species tree
        t = Tree('data/initial-data/species_tree.nw', format=1)
        orgsfile = pd.read_csv('data/initial-data/organisms.csv', header=None)

        # Create a list name/variable and use list()
        orgs = list(orgsfile[0])
        organismslist = FormatList(orgs)

        # Import alignment file as string
        alignment_file = open('data/clustal-output/' + gene + '_Aligned/' +
                              gene + '_aligned_cds_nucl.fasta', 'r')
        alignment_str = alignment_file.read()
        alignment_file.close()

        # Keep branches in the species tree for species in the alignment file
        # Some species may not be present in the alignment file
        try:
            branches2keep = []
            for organism in organismslist:
                if organism in alignment_str:
                    branches2keep.append(organism)
                else:
                    print('No sequence for %s.' % organism)

            # Input a list of branches to keep on the base tree
            t.prune(branches2keep, preserve_branch_length=True)

            # Write the tree to a file
            t.write(outfile='/work2/vallender/Projects/KARG-Project/data/paml-output/' + gene + '_PAML/temptree.nw')

            # Import the newick tree
            tree = EvolTree('/work2/vallender/Projects/KARG-Project/data/paml-output/' + gene + '_PAML/temptree.nw')

            # Import the alignment
            tree.link_to_alignment('data/clustal-output/' + gene + '_Aligned/' +
                                   gene + '_aligned_cds_nucl.fasta')

            tree.workdir = workdir

            # Set the binpath of the codeml binary
            tree.execpath = paml_path

            tree.run_model(model + '.' + gene)  # Run the model M1 M2 M3 M0

        except Exception:
            # TODO: Write a better exception.
            print('Error with %s.' % gene)
            pass
