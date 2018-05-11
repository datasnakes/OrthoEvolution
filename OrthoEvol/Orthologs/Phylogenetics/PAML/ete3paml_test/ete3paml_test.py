"""Test for codeml. Ensures codeml is working."""
from ete3 import EvolTree


class PamlTest(object):
    """Test codeml with a default tree and newick file."""

    def __init__(self, tree="ECP_EDN_15.nw", alignment="ECP_EDN_15.fasta",
                 workdir="", pamlpath=""):
        """Test that paml is in your path and working properly.

        :param tree: (Default value = "ECP_EDN_15.nw")
        :param alignment:  (Default value = "ECP_EDN_15.fasta")
        :param workdir:  (Default value = "")
        :param pamlpath:  (Default value = "")
        """

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
