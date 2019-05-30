"""Optimized for use with local/standalone NCBI BLAST 2.6.0."""
# Standard Library
import os
import shutil
import contextlib
from pathlib import Path
from datetime import datetime as d
from subprocess import run, PIPE, CalledProcessError
# BioPython
from Bio.Application import ApplicationError
from Bio import SearchIO  # Used for parsing and sorting XML files.
from Bio.Blast.Applications import NcbiblastnCommandline
# OrthoEvol
from OrthoEvol.Orthologs.Blast.base_blastn import BaseBlastN
# Other
from xml.etree.ElementTree import ParseError

# TODO-ROB: Find packages for script timing and analysis
# TODO-ROB:  Rework the save_data parameter.
# TODO-ROB:  Rework the query organism stuff.


class OrthoBlastN(BaseBlastN):
    """Combines Project Management features with NCBI's Blast+."""

    def __init__(self, project="orthology-gpcr", method=3, template=None,
                 save_data=True, acc_file="gpcr.csv", copy_from_package=True,
                 **kwargs):
        """This class inherits from the BaseBlastN class.

        This class utilizes it's parent classes to search a standalone
        Blast database for specific orthologs of a gene using a query organism
        (usually human).  The best hits from the Blast are filtered for the
        best option in order to get the most accuarate accession numbers for
        downstream analysis.

        :param project:  The project name (Default: 'orthology-inference')
        :param method: Method used for blasting. (Default: 3)
        :param template:  The accession file template.
        :param save_data:  A flag for saving the post_blast data to an excel file.
        :param acc_file: The accession file to use. (Default: 'karg.csv')
        :param copy_from_package: Copy the acc_file from the package. (Default: True)
        :param kwargs:
        """
        # Set values for methods to prevent using a config.
        self.taxon_file = None
        self.__post_blast = True
        self.go_list = None
        self.project_path = None
        self.proj_mana = None
        self.acc_file = self.MAF = acc_file
        self.copy_from_package = copy_from_package

        # Initialize class
        super().__init__(project=project, method=method, template=template,
                         save_data=save_data, acc_file=self.acc_file,
                         copy_from_package=self.copy_from_package,
                         MAF=self.MAF,
                         taxon_file=self.taxon_file,
                         post_blast=self.__post_blast,
                         go_list=self.go_list, project_path=self.project_path,
                         proj_mana=self.proj_mana, **kwargs)

    def run(self):
        """Run the blast using a default configuration."""
        self.configure(self.blast_human, self.species, auto_start=True)
