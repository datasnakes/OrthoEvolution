from Bio.Phylo.PAML import codeml
import os
from pathlib import Path


class CodemlRun(object):

    def __init__(self, P2N_alignment, iqtree_newick, control_file, home=os.getcwd()):
        self.home = Path(home)
        self.paml_path = self.home / Path('PAML')
        self.paml_path.mkdir(exist_ok=True)
        self.gene = str(iqtree_newick).replace('_iqtree.nwk', '')

        self.cml = codeml.Codeml(P2N_alignment, iqtree_newick, working_dir=home, out_file=self.gene +'_codeml.out')
        self.control_setup(control_file)

    def control_setup(self, control_file):
        self.cml.read_ctl_file(control_file)
        self.cml.print_options()
        self.cml.ctl_file = str(self.paml_path / Path(self.gene + '.ctl'))
        self.cml.write_ctl_file()



