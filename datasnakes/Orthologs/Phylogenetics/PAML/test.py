# -*- coding: utf-8 -*-
def PamlTest(pamlpath=''):
    """
    The path codeml binary should be used as the pamlpath.
    If the codeml binary is set in the environmental variables as 'codeml',
    set the pamlpath to 'codeml'.
    """
    from ete3 import EvolTree
    # Import the newick tree
    tree = EvolTree("ete3paml_test\ECP_EDN_15.nw")

    # Import the alignment
    tree.link_to_alignment("ete3paml_test\ECP_EDN_15.fasta")

    tree.workdir = 'ete3paml_test'

    # Set the binpath of the codeml binary
    tree.execpath = pamlpath

    # Run the codeml model
    tree.run_model('M1.ECP')
