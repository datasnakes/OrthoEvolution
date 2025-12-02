# -*- coding: utf-8 -*-
# Copyright 2018 by Rob Gilmore and Shaurita Hutchins. All rights reserved.
# Based on ClustalOmega wrapper copyright 2011 by Andreas Wilm.
#
# Wrapper for PAL2NAL by Rob Gilmore (2018).  http://www.bork.embl.de/pal2nal/
# Used _ClustalOmega.py as template.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Command line wrapper for PAL2NAL.
It converts a multiple sequence alignment of proteins and
the corresponding DNA (or mRNA) sequences into a codon alignment.
"""


from __future__ import print_function

import os
from Bio.Application import _Option, _Switch, AbstractCommandline, _Argument


class PAL2NALCommandline(AbstractCommandline):
    """Command line wrapper for PAL2NAL.
    http://www.bork.embl.de/pal2nal/
    Notes
    -----
    Last checked against version: v14
    References
    ----------
    Mikita Suyama, David Torrents, and Peer Bork (2006)
    PAL2NAL: robust conversion of protein sequence alignments into the
    corresponding codon alignments.
    Nucleic Acids Research, 1 July 2006; 34: W609–W612;
    https://doi.org/10.1093/nar/gkl315
    Examples
    --------
    >>> from Bio.Align.Applications import PAL2NALCommandline
    """

    def __init__(self, cmd='pal2nal', **kwargs):
        """Initialize PAL2NAL command line wrapper.

        :param cmd: Command name for PAL2NAL executable.
        :type cmd: str
        :param kwargs: Parameters for PAL2NAL configuration including:
            - pepaln: Protein alignment file (required)
            - nucfasta: DNA sequences file (required)
            - output: Output format (clustal|paml|fasta|codon)
            - output_file: Output file path (required)
            - nogap: Remove gaps and stop codons (switch)
            - nomismatch: Remove mismatched codons (switch)
            - codontable: Genetic code table number (optional)
        :type kwargs: dict
        """
        self.parameters = \
            [
                # Required Parameters
                _Argument(['pepaln'],
                          'protein alignment either in CLUSTAL or FASTA format',
                          filename=True, is_required=True,
                          checker_function=lambda x: os.path.isfile(x)),
                _Argument(['nucfasta'],
                          'DNA sequences (single multi-fasta or separated files)',
                          filename=True, is_required=True,
                          checker_function=lambda x: os.path.isfile(x)),
                _Switch(['-h', 'help'],
                        'Show help'),
                _Option(['-output', 'output'],
                        "Output format (clustal|paml|fasta|codon); default = clustal",
                        equate=False,
                        checker_function=lambda x: x in ['clustal', 'paml', 'fasta', 'codon']),
                _Switch(['-blockonly', 'blockonly'],
                        "Show only user specified blocks '#' under CLUSTAL alignment (see example)"),
                _Switch(['-nogap', 'nogap'],
                        "Remove columns with gaps and inframe stop codons"),
                _Switch(['-nomismatch', 'nomismatch'],
                        "Remove mismatched codons (mismatch between pep and cDNA) from the output"),
                _Option(['-codontable', 'codontable'],
                        "   1  Universal code (default)\
                            2  Vertebrate mitochondrial code\
                            3  Yeast mitochondrial code\
                            4  Mold, Protozoan, and Coelenterate Mitochondrial code\
                               and Mycoplasma/Spiroplasma code\
                            5  Invertebrate mitochondrial\
                            6  Ciliate, Dasycladacean and Hexamita nuclear code\
                            9  Echinoderm and Flatworm mitochondrial code\
                            10  Euplotid nuclear code\
                            11  Bacterial, archaeal and plant plastid code\
                            12  Alternative yeast nuclear code\
                            13  Ascidian mitochondrial code\
                            14  Alternative flatworm mitochondrial code\
                            15  Blepharisma nuclear code\
                            16  Chlorophycean mitochondrial code",
                        equate=False,
                        checker_function=lambda x: isinstance(x, int)),
                _Option(['>', 'output_file'],
                        "This issues the bash command that redirects the PAL2NAL"
                        "alignment to a particular file",
                        filename=True, equate=False, is_required=True)
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
