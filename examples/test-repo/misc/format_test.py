# -*- coding: utf-8 -*-
"""
Date created: Thu Mar  9 12:11:43 2017
Author: S. Hutchins

Script description: Take aligned sequences in fasta format and convert to
relaxed phylip format.

"""

from Bio import AlignIO

# Convert the file to relaxed-phylip format
AlignIO.convert("params6.fasta", "fasta", "params6.phy", "phylip-relaxed")