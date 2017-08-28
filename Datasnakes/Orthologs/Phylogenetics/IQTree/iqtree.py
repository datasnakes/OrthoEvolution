# -*- coding: utf-8 -*-
# Copyright 2017 by Rob Gilmore and Shaurita Hutchins. All rights reserved.
# Based on ClustalOmega wrapper copyright 2011 by Andreas Wilm.
#
# Wrapper for IQTree by Rob Gilmore (2017).  http://www.iqtree.org/
# Used _ClustalOmega.py as template.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Command line wrapper for IQ-Tree
Key Features:

# Efficient search algorithm: Fast and effective stochastic algorithm to reconstruct phylogenetic trees by maximum likelihood.
IQ-TREE compares favorably to RAxML and PhyML in terms of likelihood while requiring similar amount of computing time
(Nguyen et al., 2015).

# Ultrafast bootstrap: An ultrafast bootstrap approximation (UFBoot) to assess branch supports. UFBoot is 10 to 40 times
faster than RAxML rapid bootstrap and obtains less biased support values (Minh et al., 2013).

# Ultrafast model selection: An ultrafast and automatic model selection (ModelFinder) which is 10 to 100 times faster than
jModelTest and ProtTest. ModelFinder also finds best-fit partitioning scheme like PartitionFinder.

# Big Data Analysis: Supporting huge datasets with thousands of sequences or millions of alignment sites via checkpointing,
safe numerical and low memory mode. Multicore CPUs and parallel MPI system are utilized to speedup analysis.

# Phylogenetic testing: Several fast branch tests like SH-aLRT and aBayes test (Anisimova et al., 2011) and tree topology
tests like the approximately unbiased (AU) test (Shimodaira, 2002).

The strength of IQ-TREE is the availability of a wide variety of phylogenetic models:

# Common models: All common substitution models for DNA, protein, codon, binary and morphological data with rate
heterogeneity among sites and ascertainment bias correction for e.g. SNP data.

# Partition models: Allowing individual models for different genomic loci (e.g. genes or codon positions), mixed data types,
mixed rate heterogeneity types, linked or unlinked branch lengths between partitions.

# Mixture models: fully customizable mixture models and empirical protein mixture models and.
Polymorphism-aware models: Accounting for incomplete lineage sorting to infer species tree from genome-wide population
data (Schrempf et al., 2016)..
"""

from __future__ import print_function
from Bio.Application import _Option, AbstractCommandline


class IQTreeCommandline(AbstractCommandline):
    u""""Command line wrapper for GUIDANCE2.
    http://guidance.tau.ac.il/ver2/
    Example:
    --------

    \>>> from Bio.Align.Applications import IQTreeCommandline

    You would typically run the command line with clustalomega_cline() or via
        the Python subprocess module, as described in the Biopython tutorial.
    Citation:
    ---------
        To maintain IQ-TREE, support users and secure fundings, it is important for us that you cite the following papers,
        whenever the corresponding features were applied for your analysis.

        Example 1: We obtained branch supports with the ultrafast bootstrap (Minh et al., 2013) implemented in the
                    IQ-TREE software (Nguyen et al., 2015).
        Example 2: We inferred the maximum-likelihood tree using the edge-linked partition model in
                    IQ-TREE (Chernomor et al., 2016; Nguyen et al., 2015).
        ################################################################################################################
        # If you used ModelFinder please cite:
            S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, and L.S. Jermiin (2017) ModelFinder: Fast Model
            Selection for Accurate Phylogenetic Estimates, Nature Methods, 14:587–589.

        # If you performed tree reconstruction please cite:
            L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, and B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic
            algorithm for estimating maximum likelihood phylogenies. Mol. Biol. Evol., 32:268-274. DOI: 10.1093/molbev/msu300

        # If you used partition models e.g., for phylogenomic analysis please cite:
            O. Chernomor, A. von Haeseler, and B.Q. Minh (2016) Terrace aware data structure for phylogenomic inference
            from supermatrices. Syst. Biol., 65:997-1008. DOI: 10.1093/sysbio/syw037

        # If you performed the ultrafast bootstrap (UFBoot) please cite:
            B.Q. Minh, M.A.T. Nguyen, and A. von Haeseler (2013) Ultrafast approximation for phylogenetic bootstrap.
            Mol. Biol. Evol., 30:1188-1195. DOI: 10.1093/molbev/mst024

        # If you used the polymorphism-aware models please cite:
            D. Schrempf, B.Q. Minh, N. De Maio, A. von Haeseler, and C. Kosiol (2016) Reversible polymorphism-aware phylogenetic
            models and their application to tree inference. J. Theor. Biol., 407:362–370. DOI: 10.1016/j.jtbi.2016.07.042

        # If you used the IQ-TREE web server please cite:
            J. Trifinopoulos, L.-T. Nguyen, A. von Haeseler, and B.Q. Minh (2016) W-IQ-TREE: a fast online phylogenetic tool
            for maximum likelihood analysis. Nucleic Acids Res., 44 (W1):W232-W235. DOI: 10.1093/nar/gkw256
    """

    def __init__(self, cmd="iqtree", **kwargs):
        self.parameters = \
            [
                _Option(['-s', 'alignment'],
                        "Input alignment in PHYLIP/FASTA/NEXUS/CLUSTAL/MSF format",
                        filename=True, equate=False,
                        is_required=True
                        ),
                _Option(['-st', 'dataType'],
                        "BIN, DNA, AA, NT2AA, CODON, MORPH (default: auto-detect)",
                        equate=False,
                        checker_function=lambda x: x in ['BIN', 'DNA', 'AA', 'NT2AA', 'CODON', 'MORPH', 'auto-detect']),
                _Option(['', 'opts'],
                        "A placeholder to set additional parameters."
                        "e.g.  -m <model-name> -o <outgroup_taxon> -quiet -safe -mem RAM",
                        equate=False)
            ]

        AbstractCommandline.__init__(self, cmd, **kwargs)
