# -*- coding: utf-8 -*-
"""
File Name: OrthoPhyml.py
Description:

Author: S. Hutchins
Date Created: Thu May  4 16:18:08 2017
Project Name: Orthologs Project
"""

# NOTES!!! MORE DEVELOPMENT WILL BE ADDED HERE

from Bio import AlignIO
from Bio.Phylo.Applications import PhymlCommandline
import sys

class PhyML(object):
    def relaxphylip(gene):
        """Convert the file to relaxed-phylip format."""
        AlignIO.convert(gene + "_aligned_cds_nucl.fasta", "fasta",
                        gene + "_aligned.phy", "phylip-relaxed")

    def runphyml(gene):
        """Run phyml to generate tree results."""
        # Use the phyml executable file
        phyml_exe = None

        # This is mainly intended for windows use or use with an executable file
        exe_name = "PhyML-3.1_win32.exe" if sys.platform == "win32" else "phyml"
        phyml_exe = exe_name

        # Create the command & run phyml
        # Input a phylip formatted alignment file and describe the datatype ('nt' or 'aa')
        run_phyml = PhymlCommandline(phyml_exe, input=gene + '_aligned.phy', datatype='nt')
        print(run_phyml)
        out_log, err_log = run_phyml()

