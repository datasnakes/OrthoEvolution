# Standard Library
import os
from pathlib import Path
from shutil import copy
import pkg_resources
# BioPython
from Bio.Phylo.PAML import codeml
# OrthoEvol
from OrthoEvol.Manager.config import paml_control_files


class CodemlRun(object):
    """Run PAML's codeml program using Biopython's codeml wrapper.

    This class sets up the necessary files and configuration for running
    codeml analysis on a PAL2NAL alignment and IQTree newick tree.
    """

    def __init__(self, P2N_alignment, iqtree_newick, control_file='codeml-8-11-2017.ctl',
                 home=os.getcwd()):
        """Initialize CodemlRun class.

        :param P2N_alignment: Path to PAL2NAL alignment file.
        :type P2N_alignment: str
        :param iqtree_newick: Path to IQTree newick tree file.
        :type iqtree_newick: str
        :param control_file: Name of the codeml control file template.
        :type control_file: str
        :param home: Working directory for codeml execution.
        :type home: str or Path
        """
        # TODO: Generalize API and functions.

        # Set up paths
        self.home = Path(home)
        self.paml_path = self.home / Path('PAML')
        self.paml_path.mkdir(exist_ok=True)

        # Set up genes control file name and get the OrthoEvol control file path
        self.gene = str(iqtree_newick).replace('_iqtree.nwk', '')
        self.control_file = self.paml_path / Path(self.gene + '.ctl')
        self.control_template = pkg_resources.resource_filename(
            paml_control_files.__name__, control_file)
        print(self.control_template)

        # Set up CODEML input files
        self.P2N_alignment = self.home / Path(P2N_alignment)
        self.iqtree_newick = self.home / Path(iqtree_newick)
        self.P2N_alignment = copy(str(self.P2N_alignment), str(self.paml_path))
        self.iqtree_newick = copy(str(self.iqtree_newick), str(self.paml_path))
        os.chdir(str(self.paml_path))

        self.cml = codeml.Codeml(self.P2N_alignment, self.iqtree_newick,
                                 working_dir=str(self.paml_path),
                                 out_file=self.gene + '_codeml.out')
        self.control_setup(self.control_template)

    def control_setup(self, control_template):
        """Set up the codeml control file from template.

        Reads the control file template, prints options, and writes
        a gene-specific control file.

        :param control_template: Path to the codeml control file template.
        :type control_template: str
        """
        self.cml.read_ctl_file(control_template)
        self.cml.print_options()
        self.cml.ctl_file = str(self.control_file)
        self.cml.write_ctl_file()
