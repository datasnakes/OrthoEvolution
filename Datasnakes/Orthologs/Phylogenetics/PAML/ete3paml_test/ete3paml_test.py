"""Test for codeml. Ensures codeml is working."""
from ete3 import EvolTree


class PamlTest(object):
    """The path to the codeml binary should be used as the pamlpath.

    If the codeml binary is in your $PATH as 'codeml', the pamlpath should be
    the same as the default parameters.
    """
    def __init__(self, tree="ECP_EDN_15.nw", alignment="ECP_EDN_15.fasta",
                 workdir='',
                 pamlpath=''):

        self.tree = tree
        self.alignment = alignment
        self.pamlpath = pamlpath

        model = 'M1'
        self.defaultmodel = model

        wd = workdir
        self.workdir = wd

    def main(self):
        """The main function for running the test."""

        print("Running model %s paml on input." % str(self.defaultmodel))

        tree = EvolTree(self.tree)  # Import the newick tree
        tree.link_to_alignment(self.alignment)  # Import the alignment
        tree.workdir = self.workdir  # Set the working directory
        tree.execpath = self.pamlpath  # Set the binpath of the codeml binary
        tree.run_model(self.defaultmodel)  # Run the codeml model


if __name__ == "__main__":
    PamlTest().main()
