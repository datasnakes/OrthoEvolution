import os
import sys
import zipfile
import logging as log
from pathlib import Path
from datetime import datetime as d
from Datasnakes.Manager.utils.mana import ProjMana
from Datasnakes.Orthologs.Blast.blastn import BLASTn
from Datasnakes.Orthologs.CompGenetics.comp_gen import CompGenAnalysis
from Datasnakes.Orthologs.Genbank.genbank import GenBank
from Datasnakes.Orthologs.Align.msa import MultipleSequenceAlignment as MSA
import yaml
#import configparser
#from slacker import Slacker
#import argparse
#import textwrap
# TODO-ROB: Use this class to move data back and forth between local and remote servers
# TODO-ROB:  Mirror the directory creation on MCSR's servers
# TODO-ROB:  ^^ This will allow the transfer of data

# TODO-ROB:  Add FTP and s2s inheritance
class DataMana(object):

    def __init__(self, config_file=None, pipeline=None, **kwargs):

        self.ProjectManagment_config = self.CompGenAnalysis_config = self.BLASTn_config = self.GenBank_config = self.Alignment_config = None

        if pipeline == 'Ortho_CDS_1':
            config_file = pkg_resources.resource_filename()

        if config_file is not None:
            with open(config_file, 'r') as ymlfile:
                configuration = yaml.load(ymlfile)
                for key, value in configuration.items():
                    setattr(self, key, value)

            # Project Management
            if self.ProjectManagment_config is not None:
                self.pm = ProjMana(**self.ProjectManagment_config)
            else:
                self.pm = self.ProjectManagment_config

            # BLASTn
            if self.BLASTn_config is not None:
                self.blast(self.pm)
            else:
                self.bl = self.BLASTn_config

            # GenBank
            if self.GenBank_config is not None:
                self.genbank(self.bl)
            else:
                self.gb = self.GenBank_config

            # Alignment
            if self.Alignment_config is not None:
                self.align(self.gb)
            else:
                self.al = self.Alignment_config

    def blast(self, proj_mana):
        self.bl = BLASTn(proj_mana=proj_mana, **self.BLASTn_config)
        self.bl.blast_config(self.bl.blast_human, 'Homo_sapiens', auto_start=True)
        # TODO-Create directories for the blast data
        # Do the blasting here using BLASTn

    def genbank(self, blast):
        self.gb = GenBank(blast=blast, **self.GenBank_config)
        if isinstance(blast, BLASTn):
            self.gb.blast2_gbk_files(blast.org_list, blast.gene_dict)
        else:

            cga = CompGenAnalysis(**self.CompGenAnalysis_config)

            # Parse the tier_frame_dict to get the tier
            for G_KEY, G_value in cga.tier_frame_dict.items():
                tier = G_KEY
                # Parse the tier based transformed dataframe to get the gene
                for GENE in cga.tier_frame_dict[tier].T:
                    # Parse the organism list to get the desired accession number
                    for ORGANISM in cga.org_list:
                        accession = str(cga.gene_dict[GENE][ORGANISM])
                        accession, sup, version = accession.partition('.')
                        accession = accession.upper()
                        server_flag = False
                        self.bl.get_gbk_file(accession, GENE, ORGANISM, server_flag=server_flag)

    def align(self, genbank):
        self.al = MSA(genbank=genbank, **self.Alignment_config)
        self.al.align(self.Alignment_config['kwargs'])


class ZipUtils(DataMana):
    """The ZipUtil class allows easy compression/zipping of file folders.
    Inspired by http://stackoverflow.com/a/670635/7351746
    """

    def __init__(self, zip_filename, zip_path):
        """Initialize the input files and path.

        :param comp_filename (string):  This is the name of the compressed file that will be generated (eg 'test.zip')
        :param zip_path: This is the absolute path of the directory (or file) to be zipped.
        :returns:  A zip file that is created inside of the zip_path.  The path string is returned.
        """
        self.zip_filename = zip_filename
        self.zip_path = zip_path
        self.ignore_parts = Path(zip_path).parent.parts

    def to_zip(self):
        """Zip a folder."""
        comp_path = os.path.join(self.zip_path, self.zip_filename)
        zip_handle = zipfile.ZipFile(comp_path, 'w', zipfile.ZIP_DEFLATED)
        if os.path.isfile(self.zip_path):
            zip_handle.write(self.zip_path)
        else:
            print('skipped')
            self.add_folder_to_zip(zip_handle, self.zip_path)
        zip_handle.close()
        return comp_path

    def add_folder_to_zip(self, zip_handle, folder):
        """Not meant to be used explicitly.  Use to_zip.

        :param zip_handle: An initialized zipfile.ZipFile handle.
        :param folder: A path that represents an entire folder to be zipped.
        :return: Recursively zips nested directories.
        """
        for file in os.listdir(folder):
            full_path = os.path.join(folder, file)
            rel_path = Path(full_path)
            rel_path = rel_path.relative_to(Path(self.zip_path))
            if os.path.isfile(full_path):
                if str(file) == str(self.comp_filename):
                    continue
                print('File added: ' + str(full_path))
                zip_handle.write(full_path, rel_path)
            elif os.path.isdir(full_path):
                if str(file) in self.ignore_parts:
                    continue
                print('Entering folder: ' + str(full_path))
                self.add_folder_to_zip(zip_handle, full_path)
