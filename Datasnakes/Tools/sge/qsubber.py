import argparse
import textwrap
from qstat import Qstat

__author__ = 'Datasnakes'

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
                                    This is a command line wrapper for the SGE module.

                                    ' '''))
parser.add_argument("-o", "--output", help="Qstat info output type",
                    required=True)

q = Qstat()
args = parser.parse_args(namespace=q)
