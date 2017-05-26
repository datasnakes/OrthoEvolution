# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 15:22:52 2017

@author: S. Hutchins
"""

# Translate nucleotide cds to amino acid sequence

# Modules used
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Import coding dna
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", generic_dna)

# Translate the dna & save to a separate file
with open("gene_cds_aa.txt", "w") as aa:
    aa.write(str(">" + "\n" + coding_dna.translate(table=1, to_stop=True)))





