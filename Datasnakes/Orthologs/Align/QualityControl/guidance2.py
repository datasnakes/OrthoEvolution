# -*- coding: utf-8 -*-
# Copyright 2017 by Rob Gilmore. All rights reserved.
# Based on ClustalOmega wrapper copyright 2011 by Andreas Wilm.
#
# Wrapper for Guidance2 by Rob Gilmore (2017).  http://guidance.tau.ac.il/ver2/
# Used _ClustalOmega.py as template.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the weighting, filtering or masking
of unreliably aligned positions in multiple sequence alignments using GUIDANCE2."""

from __future__ import print_function
from pathlib import Path
from Bio.Application import _Option, _Switch, AbstractCommandline

class Guidance2Commandline(AbstractCommandline):
    u""""Command line wrapper for GUIDANCE2.
    http://guidance.tau.ac.il/ver2/
    Example:
    --------

    \>>> from Bio.Align.Applications import Guidance2Commandline

    You would typically run the command line with clustalomega_cline() or via
        the Python subprocess module, as described in the Biopython tutorial.
    Citation:
    ---------
    Sela, I., Ashkenazy, H., Katoh, K. and Pupko, T. (2015)
    GUIDANCE2: accurate detection of unreliable alignment regions accounting for the uncertainty of multiple parameters.
    Nucleic Acids Research, 2015 Jul 1; 43 (Web Server issue): W7-W14.; doi: 10.1093/nar/gkq443

    Landan, G., and D. Graur. (2008).
    Local reliability measures from sets of co-optimal multiple sequence alignments.
    Pac Symp Biocomput 13:15-24
    """

    def __init__(self, cmd="guidance", **kwargs):
        # order parameters in the same order as invoking guidance on the cmd line (e.g. 'perl guidance.pl')
        self.parameters = \
            [
                # Required Parameters
                _Option(['--seqFile', 'seqFile'],
                        "Input sequence file in FASTA format",
                        filename=True, equate=False, is_required=True,
                        checker_function=lambda x: Path(x).suffix in ['.fasta', 'fna', '.ffn', '.faa', '.fra']),
                _Option(['--msaProgram', 'msaProgram'],
                        "Which MSA program to use",
                        equate=False, is_required=True,
                        checker_function=lambda x: x in ['MAFFT', 'PRANK', 'CLUSTALW', 'MUSCLE']),
                _Option(['--seqType', 'seqType'],
                        "Type of sequences for alignment (amino acids, nucleotides, or codons)",
                        equate=False, is_required=True,
                        checker_function=lambda x: x in ['aa', 'nuc', 'codon']),
                _Option(['--outDir', 'outDir'],
                        "Output directory that will be created "
                        "automatically and hold all output files [please provid full (and not relative) path]",
                        filename=True, equate=False, is_required=True,
                        checker_function=lambda x: Path(x).is_dir()),

                # Optional Parameters
                _Option(['--program', 'program'],
                        "[GUIDANCE2|GUIDANCE|HoT] Default=GUIDANCE2",
                        equate=False,
                        checker_function=lambda x: x in ["GUIDANCE2", "GUIDANCE", "HoT"]),
                _Option(['--bootstraps', 'bootstraps'],
                        "",
                        equate=False),
                _Option(['--genCode', 'genCode'],
                        "Genetic code identifier (only for codon sequences). Default=1 \
                            1) Nuclear Standard\
                            15) Nuclear Blepharisma\
                            6) Nuclear Ciliate\
                            10) Nuclear Euplotid\
                            2) Mitochondria Vertebrate\
                            5) Mitochondria Invertebrate\
                            3) Mitochondria Yeast\
                            13) Mitochondria Ascidian\
                            9) Mitochondria Echinoderm\
                            14) Mitochondria Flatworm\
                            4) Mitochondria Protozoan",
                        equate=False),
                _Option(['--outOrder', 'outOrder'],
                        "",
                        equate=False),
                _Option(['--msaFile', 'msaFile'],
                        "",
                        equate=False),
                _Option(['', ''],
                        "",
                        equate=False),
                _Option(['', ''],
                        "",
                        equate=False),
                _Option(['', ''],
                        "",
                        equate=False),
                _Option(['', ''],
                        "",
                        equate=False),
                _Option(['', ''],
                        "",
                        equate=False),
                _Option(['', ''],
                        "",
                        equate=False),
                _Option(['', ''],
                        "",
                        equate=False),


            ]
