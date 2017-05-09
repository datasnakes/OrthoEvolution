# -*- coding: utf-8 -*-
from ete3 import EvolTree, Tree
import pandas as pd

## Import the newick tree
#tree = EvolTree("Ete3Paml-Examples\species_tree.txt")
#
## Import the alignment
#tree.link_to_alignment("Ete3Paml-Examples\HTR1D.fas")
#
#tree.workdir = 'Ete3Paml-Examples'
#
## Set the binpath of the codeml binary
#tree.execpath = r'C:\Users\shutchins2\Desktop\Software & Executables\paml4.9c\bin'
#
## Run the codeml model
#tree.run_model('M1.HTR1D')
#def prune_tree():
t = Tree('Ete3Paml-Examples\species_tree.txt', format=1)

# Create a list of
#orgsfile = pd.read_csv('data/initial-data/organisms.csv', header=None)
orgsfile = pd.read_csv('Ete3Paml-Examples\organisms.csv', header=None)

# Create a list name/variable and use list()
orgs = list(orgsfile[0])
def formatlist(input_list):
    """Remove spaces from list items and turn those spaces into underscores."""
    output_list = []
    for item in input_list:
        item = str(item)
        item = item.replace(" ", "_")
        output_list.append(item)
        return output_list
organismslist = formatlist(orgs)

# Import alignment file as string
alignment_file = open('Ete3Paml-Examples\HTR1D.fas', 'r')
alignment_str = alignment_file.read()

branches2keep = []
for organism in organismslist:
    if organism in alignment_str:
        print('Yup.')
        branches2keep.append(organism)
    else:
        print('Nope.')


# Input a list of branches to keep on the base tree
t.prune(branches2keep, preserve_branch_length=True)