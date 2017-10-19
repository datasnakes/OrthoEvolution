import argparse
import textwrap
from qstat import Qstat

__author__ = 'Datasnakes'

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
                                    This is a command line wrapper for the sge module.

                                    ' '''))
parser.add_argument("-o", "--output", help="Qstat info output type",
                    required=False)

q = Qstat()
args = parser.parse_args(namespace=q)

# TODO Complete this