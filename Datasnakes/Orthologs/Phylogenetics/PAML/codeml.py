from Bio.Phylo.PAML import codeml
import os
from pathlib import Path
from shutil import copy
from Datasnakes.Orthologs.Phylogenetics.PAML import ctlfiles
control_path = ctlfiles.__path__._path[0]


class CodemlRun(object):

    def __init__(self, P2N_alignment, iqtree_newick, control_file='codeml-8-11-2017.ctl', home=os.getcwd()):
        # Set up paths
        self.home = Path(home)
        self.paml_path = self.home / Path('PAML')
        self.paml_path.mkdir(exist_ok=True)
        # Set up genes control file name and get the datasnakes control file path
        self.gene = str(iqtree_newick).replace('_iqtree.nwk', '')
        self.control_file = self.paml_path / Path(self.gene + '.ctl')
        self.control_template = Path(control_path) / Path(control_file)

        # Set up CODEML input files
        self.P2N_alignment = self.home / Path(P2N_alignment)
        self.iqtree_newick = self.home / Path(iqtree_newick)
        self.P2N_alignment = copy(str(self.P2N_alignment), str(self.paml_path))
        self.iqtree_newick = copy(str(self.iqtree_newick), str(self.paml_path))
        os.chdir(str(self.paml_path))

        self.cml = codeml.Codeml(self.P2N_alignment, self.iqtree_newick, working_dir=str(self.paml_path), out_file=self.gene +'_codeml.out')
        self.control_setup(self.control_template)

    def control_setup(self, control_template):
        self.cml.read_ctl_file(control_template)
        self.cml.print_options()
        self.cml.ctl_file = str(self.control_file)
        self.cml.write_ctl_file()



