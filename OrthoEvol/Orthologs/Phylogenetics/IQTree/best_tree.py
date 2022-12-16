# Standard Library
import os
from shutil import copy
from subprocess import check_call, STDOUT
from pathlib import Path
# OrthoEvol
from OrthoEvol.Orthologs.Phylogenetics.IQTree.iqtree import IQTreeCommandline
#TODO-ROB Make this inherit FilteredAlignment

# Set up make directory function


class FilteredTree(object):
    """Run IQTree to generate a "filtered" tree or best tree.

    :param alignment: Path to multiple sequence alignment file.
    :param dataType:  Input datatype. (Default value = 'CODON')
    :param home: Path of working directory.  (Default value = PWD)
        """

    def __init__(self, alignment, dataType='CODON', home=os.getcwd()):
        self.home = Path(home)
        self.iqtree_path = self.home / Path('IQTREE')
        self.tree_file = self.iqtree_path / Path(alignment + '.treefile')
        self.gene = alignment.replace('_P2N_na.iqtree.aln', '')

        self.aln_File = str(self.home / Path(alignment))
        outDir = self.home / Path('IQTREE')
        os.makedirs(str(outDir), exist_ok=True)
        copy(self.aln_File, str(outDir))
        os.chdir(str(outDir))
        self.iqtree_best_tree(alignment, dataType)
        copy(self.tree_file, str(self.home / Path(self.gene + '_iqtree.nwk')))

    def iqtree_best_tree(self, alignment, dataType):
        """Generate and save the best tree from IQTree by importing the filtered
        alignment.
        """
        iqtree_cline = IQTreeCommandline(alignment=alignment,
                                         dataType=dataType)
        print(iqtree_cline)
        check_call([str(iqtree_cline)], stderr=STDOUT, shell=True)
