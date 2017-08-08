"""Test for codeml. Ensures codeml is working."""
from ete3 import EvolTree


class PamlTest:
    """The path to the codeml binary should be used as the pamlpath.
    
    If the codeml binary is set in the environmental variables as 'codeml',
    set the pamlpath to 'codeml'.
    """
    def __init__(pamlpath):
        # Import the newick tree
        tree = EvolTree("ECP_EDN_15.nw")
    
        # Import the alignment
        tree.link_to_alignment("ECP_EDN_15.fasta")
    
        tree.workdir = 'test'
    
        # Set the binpath of the codeml binary
        tree.execpath = pamlpath
    
        # Run the codeml model
        tree.run_model('M1.ECP')