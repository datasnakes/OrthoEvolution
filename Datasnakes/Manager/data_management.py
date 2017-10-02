import os
import zipfile
from pathlib import Path

import pkg_resources
import yaml

from Datasnakes.Manager import config
from Datasnakes.Manager import ProjectManagement
from Datasnakes.Orthologs.Align import MultipleSequenceAlignment as MSA
from Datasnakes.Orthologs.Blast.blastn_comparative_genetics import CompGenBLASTn
from Datasnakes.Orthologs.Blast.comparative_genetics_objects import CompGenObjects
from Datasnakes.Orthologs.Genbank.genbank import GenBank


#import configparser
#from slacker import Slacker
#import argparse
#import textwrap
# TODO-ROB: Use this class to move data back and forth between local and remote servers
# TODO-ROB:  Mirror the directory creation on MCSR's servers
# TODO-ROB:  ^^ This will allow the transfer of data

# TODO-ROB:  Add FTP and s2s inheritance


class DataMana(object):

    def __init__(self, config_file=None, pipeline=None, new=False, start=False, **kwargs):
        """Initialize the attributes that can be used as keys in the config_file."""
        self.Management_config = self.CompGenAnalysis_config = self.BLASTn_config = self.GenBank_config = self.Alignment_config = None
        self.pm = self.bl = self.gb = self.al = None
        if pipeline == 'Ortho_CDS_1':
            if new is True:
                config_file = pkg_resources.resource_filename(config.__name__, 'config_template_new.yml')
            else:
                config_file = pkg_resources.resource_filename(config.__name__, 'config_template_existing.yml')
        if config_file is not None:
            if start is True:
                self.configure(config_file)

    def configure(self, config_file):
        '''
        This method uses YAML configuration in order to initialize different classes.
        :param config_file: A YAML file that is used to create a dictionary(kwargs) for each class.
        :return:
        '''
        with open(config_file, 'r') as ymlfile:
            configuration = yaml.safe_load(ymlfile)
            for key, value in configuration.items():
                setattr(self, key, value)
                print('key:' + str(key) + '\nvalue: ' + str(value))

                # Project Management
            if self.Management_config is not None:
                self.pm = ProjectManagement(**self.Management_config)
                print('mana_config')
                print(self.pm)
            else:
                self.pm = self.Management_config

                # CompGenAnalysis and BLASTn Configuration
            if self.BLASTn_config is not None and self.CompGenAnalysis_config is not None:
                self.BLASTn_config.update(self.CompGenAnalysis_config)
            if self.BLASTn_config is not None:
                # Blast has not taken place so it will happen here
                self.blast(self.pm, self.BLASTn_config)
            elif self.CompGenAnalysis_config is not None:
                self.blast(self.pm, self.CompGenAnalysis_config)
            else:
                # Blast has taken place so
                self.bl = self.BLASTn_config

                # GenBank
            if self.GenBank_config is not None:
                self.genbank(self.pm, self.bl)
            else:
                self.gb = self.GenBank_config

                # Alignment
            if self.Alignment_config is not None:
                self.align(self.gb)
            else:
                self.al = self.Alignment_config

    def blast(self, proj_mana, blast_config):
        self.bl = CompGenBLASTn(proj_mana=proj_mana, **self.Management_config, **blast_config)
        self.bl.blast_config(self.bl.blast_human, 'Homo_sapiens', auto_start=True)
        # TODO-Create directories for the blast data
        # Do the blasting here using CompGenBLASTn

    def genbank(self, proj_mana, blast):
        if blast is not None:
            self.gb = GenBank(blast=blast, **self.Management_config, **self.GenBank_config)
        else:
            self.gb = GenBank(blast=blast, **self.Management_config, **self.GenBank_config)
        if blast is not None:
            if issubclass(type(blast), CompGenBLASTn):
                self.gb.blast2_gbk_files(blast.org_list, blast.gene_dict)
        else:
            print(proj_mana.__dict__)
            cga = CompGenObjects(proj_mana=proj_mana, **self.CompGenAnalysis_config)

            # Parse the tier_frame_dict to get the tier
            for G_KEY in cga.tier_frame_dict.keys():
                tier = G_KEY
                # Parse the tier based transformed dataframe to get the gene
                for GENE in cga.tier_frame_dict[tier].T:
                    # Parse the organism list to get the desired accession number
                    for ORGANISM in cga.org_list:
                        accession = str(cga.gene_dict[GENE][ORGANISM])
                        parts = list(accession.partition('.'))
                        accession = parts[0]
                        accession = accession.upper()
                        server_flag = False
                        self.gb.get_gbk_file(accession, GENE, ORGANISM, server_flag=server_flag)

    def align(self, genbank):
        self.al = MSA(genbank=genbank, **self.Management_config, **self.Alignment_config)
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
