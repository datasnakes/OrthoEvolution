import argparse
import textwrap
from Datasnakes.Tools.sge import SGEJob

__author__ = 'Datasnakes'

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
                                    This is a command line wrapper for the sge module.

                                    ' '''))
