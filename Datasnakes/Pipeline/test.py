"""Simple test of one gene through the entire pipeline."""
# import os
# import sys
# import zipfile
import luigi
import logzero as log
from pathlib import Path
from datetime import datetime as d

from Datasnakes.Manager.utils.mana import ProjMana  # Project Management
from Datasnakes.Orthologs.Blast.blastn import BLASTn
from Datasnakes.Orthologs.Genbank.genbank import GenBank
from Datasnakes.Orthologs.Align.alignment import Alignment


class Blast2PAML(luigi.Task):
    """Blasn to PAML."""
    n = luigi.Parameter(default=10)

    def requires(self):
        return [ProjMana(),
                BLASTn(),
                GenBank(),
                Alignment()]

    def output(self):
        return luigi.LocalTarget("numbers_up_to_{}.txt".format(self.n))

    def run(self):
        with self.output().open('w') as f:
            for i in range(1, self.n+1):
                f.write("{}\n".format(i))


if __name__ == '__main__':
    # from luigi.mock import MockFile
    luigi.run(["--local-scheduler"], main_task_cls='Blast2PAML')
