#!/usr/bin/env python
"""This script is designed to generate some basic gene information from a csv
file of refseqrna accession numbers for human genes."""

import argparse
import textwrap

from OrthoEvol.Tools.mygene import MyGene


def main(infile, outfile):
    """Use MyGene to generate basic gene information.

    :param infile: Path to csv input file.
    :param outfile: Path to csv output file.
    """
    mg = MyGene(infile, outfile)
    mg.query_mygene()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                        This is a command line interface the uses mygene 
                                        and pandas to convert a csv column of refseqrna 
                                        accession numbers of human genes to gene information.
                                        '''))
    parser.add_argument('-i', '--infile', help='Name and path of your input file.',
                        required=True)
    parser.add_argument('-o', '--outfile',
                        help='Name and path of your output file.',
                        required=True)

    args = parser.parse_args()

    main(args.infile, args.outfile)
