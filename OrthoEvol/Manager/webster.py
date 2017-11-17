import re
import textwrap as twrap
from pathlib import Path

from OrthoEvol.Manager.config import references
from pkg_resources import resource_filename


# Based off of https://github.com/etetoolkit/ete/blob/master/ete3/citation.py
class Webster(object):

    archive_options = {
        "Full": Path(''),
        "NCBI": Path('NCBI'),
        "ITIS": Path('ITIS'),
        "NCBI_blast": Path('NCBI/blast'),
        "NCBI_blast_db": Path('NCBI/blast/db'),
        "NCBI_blast_windowmasker_files": Path('NCBI/blast/windowmasker_files'),
        "NCBI_pub_taxonomy": Path('NCBI/pub/taxonomy'),
        "NCBI_refseq_release": Path('NCBI/refseq/release'),
        "ITIS_taxonomy": Path('ITIS/taxonomy'),
    }

    byte_size_options = {
        "B": 1,
        "KB": 1024,
        "MB": 1048576,
        "GB": 1073741824,
        "TB": 1099511627776
    }

    reference_options = {
        'GUIDANCE2': {

            "reference_1": {
                "citation": u"""Sela, I., Ashkenazy, H., Katoh, K. and Pupko, T. (2015) GUIDANCE2: Accurate Detection of 
                Unreliable Alignment Regions Accounting for the Uncertainty of Multiple Parameters.  Nucleic Acids 
                Research, 2015 Jul 1; 43 (Web Server issue): W7-W14.; doi: 10.1093/nar/gkq443""",
                "link": "https://www.ncbi.nlm.nih.gov/pubmed/18229673",
                "path": resource_filename(
                    references.__name__, "GUIDANCE2_Accurate_Detection_of_Unreliable_Alignment_Regions_Accounting_for_"
                                         "the_Uncertainty_of_Multiple_Parameters.pdf")
            },
            "reference_2": {
                "citation": u"""Landan, G., and D. Graur. (2008).  Local Reliability Measures from Sets of Co-optimal 
                Multiple Sequence Alignments.  Pac Symp Biocomput 13:15-24""",
                "link": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489236/",
                "path": resource_filename(
                    references.__name__, "Local_Reliability_Measures_from_Sets_of_Co-optimal_Multiple_Sequence_"
                                         "Alignments.pdf")
            }
        },
        'PAL2NAL': u"""Mikita Suyama, David Torrents, and Peer Bork (2006) PAL2NAL: Robust Conversion of Protein 
        Sequence Alignments Into the Corresponding Codon Alignments.  Nucleic Acids Res. 34, W609-W612.""",

        'IQTREE': u"""LT Nguyen, H.A. Schmidt, A. von Haeseler, and BQ Minh (2015) IQ-TREE: A
        fast and effective stochastic algorithm for estimating maximum likelihood
        phylogenies. Mol. Biol. Evol., 32, 268-274.""",

        'CLUSTALO': u""" Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, 
        Söding J, Thompson JD, Higgins DG.  Fast, Scalable Generation of High-quality Protein Multiple Sequence 
        Alignments Using Clustal Omega.  Mol Syst Biol. 2011 Oct 11;7:539. doi: 10.1038/msb.2011.75.""",


    }

    def __init__(self):
        self.citations = set()

    def add(self, ref):
        self.citations.add(self.reference_options[ref])

    def show(self):
        wrapper = twrap.TextWrapper(width=75, initial_indent="   ",
                                    subsequent_indent="      ",
                                    replace_whitespace=False)
        citations = sorted(self.citations)
        print("   ========================================================================")
        print("         The following published software and/or methods were used.        ")
        print("               *** Please, do not forget to cite them! ***                 ")
        print("   ========================================================================")
        for ref in citations:
            print(wrapper.fill(re.sub('[\n \t]+', ' ', ref).strip()))
