import os
from shutil import copy
from subprocess import check_call, STDOUT
from pathlib import Path

from OrthoEvol.Tools import LogIt
from OrthoEvol.Orthologs.Phylogenetics.IQTree.iqtree import IQTreeCommandline
from OrthoEvol.Tools.otherutils import makedirectory
# TODO-ROB Make this inherit FilteredAlignment


class FilteredTree(object):
    """This is a  wrapper around the IQTree wrapper to get the best tree."""

    def __init__(self, alignment, dataType='CODON', working_dir=''):
        """Run IQTree to generate a "filtered" tree or best tree.

        :param alignment: Path to multiple sequence alignment file.
        :param dataType:  Input datatype. (Default value = 'CODON')
        :param working_dir: Path of working directory.  (Default value = '')
        """
        self.iqtree_log = LogIt().default(logname="iqtree", logfile=None)
        self.working_dir = Path(working_dir)
        self.iqtree_path = self.working_dir / Path('IQTREE')
        self.tree_file = self.iqtree_path / Path(alignment + '.treefile')
        self.gene = alignment.replace('_P2N_na.iqtree.aln', '')

        self.aln_File = str(self.working_dir / Path(alignment))
        outDir = self.working_dir / Path('IQTREE')
        makedirectory(outDir)
        copy(self.aln_File, str(outDir))
        os.chdir(str(outDir))
        self.iqtree_best_tree(alignment, dataType)
        treepath = str(self.working_dir / Path(self.gene + '_iqtree.nwk'))
        copy(self.tree_file, treepath)

    def iqtree_best_tree(self, alignment, dataType):
        """Generate and save the best tree from IQTree.

        :param alignment:  Path to multiple sequence alignment file.
        :param dataType:
        :return:
        """

        iqtree_cline = IQTreeCommandline(alignment=alignment,
                                         dataType=dataType)
        self.iqtree_log.info(iqtree_cline)
        check_call([str(iqtree_cline)], stderr=STDOUT, shell=True)
