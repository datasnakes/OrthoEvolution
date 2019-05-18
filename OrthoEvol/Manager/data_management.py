# Standard Library
import pkg_resources
import yaml
# OrthoEvol
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.Manager import config
from OrthoEvol.Manager.database_management import BaseDatabaseManagement
from OrthoEvol.Orthologs.Align import MultipleSequenceAlignment as MSA
from OrthoEvol.Orthologs.Blast.orthologs_blastn import OrthoBlastN
from OrthoEvol.Orthologs.Blast.comparative_genetics import BaseComparativeGenetics
from OrthoEvol.Orthologs.GenBank.genbank import GenBank


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
        # Full configuration for the pipeline's YAML based variables
        self.Management_config = self.Database_config = self.CompGenAnalysis_config = self.BLASTn_config = \
            self.GenBank_config = self.Alignment_config = None
        # Alignment configuration
        self.Guidance_config = self.Clustalo_config = self.Pal2Nal_config = None
        self.pm = self.bl = self.gb = self.al = self.db = None
        if pipeline == 'Ortho_CDS_1':
            if new is True:
                config_file = pkg_resources.resource_filename(config.yaml.__name__, 'pipeline.yml')
            else:
                config_file = pkg_resources.resource_filename(config.yaml.__name__, 'config_template_existing.yml')
        if config_file is not None:
            if start is True:
                self.configure(config_file)

    def configure(self, config_file):
        """Use YAML configuration in order to initialize different classes.
        
        :param config_file: A YAML file that is used to create a dictionary(kwargs) for each class.
        :return:
        """
        with open(config_file, 'r') as ymlfile:
            configuration = yaml.safe_load(ymlfile)
            # TODO-ROB:  Set up configuratioin to parse this full list.  For each sub configuration do something.
            # TODO-ROB:  Before doing the above set up Airflow.
            setattr(self, "CONFIGURATION", configuration)
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

            if self.Database_config is not None:
                self.database(self.pm, self.Database_config)
            else:
                self.db = self.Database_config

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

    def database(self, proj_mana, database_config):
        self.db = BaseDatabaseManagement(proj_mana=proj_mana, **database_config)
        for config_type, database_config_list in self.db.database_dict.items():
            implementation = database_config_list[0]
            configuration = database_config_list[1]
            if isinstance(implementation, dict):
                for sub_config_type, sub_implementation in implementation.items():
                    if isinstance(configuration, dict):
                        sub_implementation(**configuration[sub_config_type])
                    else:
                        sub_implementation(**configuration)
            else:
                implementation(**configuration)


            # TODO-ROB parse the config options

    def blast(self, proj_mana, blast_config):
        self.bl = OrthoBlastN(proj_mana=proj_mana, **self.Management_config, **blast_config)
        self.bl.blast_config(self.bl.blast_human, 'Homo_sapiens', auto_start=True)
        # TODO-Create directories for the blast data
        # Do the blasting here using CompGenBLASTn

    def genbank(self, proj_mana, blast):
        if blast is not None:
            self.gb = GenBank(blast=blast, **self.Management_config, **self.GenBank_config)
        else:
            self.gb = GenBank(blast=blast, **self.Management_config, **self.GenBank_config)
        if blast is not None:
            if issubclass(type(blast), OrthoBlastN):
                self.gb.create_post_blast_gbk_records(blast.org_list, blast.gene_dict)
        else:
            print(proj_mana.__dict__)
            cga = BaseComparativeGenetics(proj_mana=proj_mana, **self.CompGenAnalysis_config)

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
        for _program, aligner_config_list in self.al.alignment_dict.items():
            aligner = aligner_config_list[0]
            configuration = aligner_config_list[1]
            aligner(**configuration)
