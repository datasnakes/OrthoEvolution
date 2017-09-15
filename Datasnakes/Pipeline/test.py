"""Simple test of one gene through the entire pipeline."""
# import os
# import sys
# import zipfile
import luigi
import logzero as log
from pathlib import Path
from datetime import datetime as d
import yaml

from Datasnakes.Manager.utils import ProjMana  # Project Management
from Datasnakes.Orthologs.Blast import BLASTn
from Datasnakes.Orthologs.Genbank import GenBank
from Datasnakes.Orthologs.Align import MultipleSequenceAlignment as MSA

#logger = log.logger('luigi-interface')

# TODO Check out https://stackoverflow.com/questions/39996544/handling-a-lot-of-parameters-in-luigi

class Blast2PAML(luigi.Task):
    """Blasn to PAML."""
    # TODO Create a task family.
    # Set luigi Paramaters. Allows command line usage.
    repository = luigi.Parameter(default='test-repo')
    username = luigi.Parameter(default='test-user')
    project = luigi.Parameter(default='test-project')
    research = luigi.Parameter(default='test-research')
    researchtype = luigi.Parameter(default='comparative_genetics')
    newrepo = luigi.Parameter(default=True)
    newuser= luigi.Parameter(default=True)
    newproject = luigi.Parameter(default=True)
    newresearch = luigi.Parameter(default=True)


    def requires(self):
        return ProjMana(repo=self.repository, user=self.username,
                                     project=self.project, research_type=self.researchtype,
                                     new_repo=self.newrepo, new_user=self.newuser,
                                     new_project=self.newproject, new_research=self.newresearch)

    def output(self):
        return luigi.LocalTarget('~/')

#    def run(self):
#        with self.input()[0]


if __name__ == '__main__':
    #from luigi.mock import MockFile
    luigi.run(["--local-scheduler"], main_task_cls='Blast2PAML')
