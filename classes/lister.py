##############################################################################
# PyCharm Community Edition
# -*- coding: utf-8 -*-
"""
GPCR-Orthologs-Project
Accession2 updated on 11/17/2016 at 1:09 PM
##############################################################################

    Input:  An open .csv file object that contains a header of organisms.  The
    first column ranks the gene by tier, the second column is a HUGO Gene
    Nomenclature Committee(HGNC) symbol for the genes of interest.  The .csv
    has to be located in the same directory as this module unless a full path is
    specified.

    The organisms are taken from
    ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/multiprocessing/
    and the genes are taken from http://www.guidetopharmacology.org/targets.jsp.

    Output:  A pandas Data-Frame, Pivot-Table, and associated lists and dictionaries.

    Description:  Parses an accession file with the designated format in order to
    provide easy handles for different pieces of data.

    Notes:  It doesn't matter what tier.  Just parse the file.

##############################################################################
@author: rgilmore
"""
##############################################################################
# Libraries:

import os
from pathlib import Path

import pandas as pd

from dir_mana import dir_mana

##############################################################################
# Directory Initializations:
# Use dir_mana() class here so that we can stay organized
# and more easily access the proper directories on command
os.chdir(os.path.dirname(__file__))
home = os.getcwd()
project = "Orthologs-Project"
where = dir_mana(home, project)

# Add a path that contains custom libraries for import
# os.sys.path.append()
##############################################################################
# Initializations:

##############################################################################


class Lister(object):
    __home = ''
    __filename = ''
    __filename_path = ''
    __data = ''
    ##########################################################################

    def __init__(self, csv_file, go_list=None, hgnc=False):

        # Private Variables
        self.__home = home
        self.__filename = csv_file
        self.__filename_path = 'C:\\Users\\rgilmore\\PycharmProjects\\Orthologs-Project\\' / Path(self.__filename)

        self.go_list = go_list
        # Handles for organism lists #
        self.org_list = []
        self.org_count = 0
        # Handles for gene lists #
        self.gene_list = []
        self.gene_count = 0
        # Handles for tier lists #
        self.tier_list = []
        self.tier_dict = {}
        self.tier_frame_dict = {}
        # Handles for accession lists #
        self.acc_dict = {}
        self.acc_list = []
        # Handles for master data structures #
        self.__data = pd.read_csv(self.__filename_path)
        self.header = self.__data.axes[1].tolist()
        # #### Format the main data frame #### #
        self.__data = self.__data.set_index('Gene')
        self.df = self.__data
        # #### Format the main pivot table #### #
        self.pt = pd.pivot_table(pd.read_csv(self.__filename_path), index=['Tier', 'Gene'], aggfunc='first')
        array = self.pt.axes[1].tolist()
        self.pt.columns = pd.Index(array, name='Organism')
        # #### Handles for full dictionaries #### #
        self.org_dict = self.df.ix[0:, 'Homo_sapiens':].to_dict()
        self.gene_dict = self.df.T.to_dict()
        self.get_master_lists(self.__data)  # populates our lists

        # TODO-ROB Add HGNC python module

    def get_accession(self, gene, organism):
        """Takes a single gene and organism and returns
        a single accession number."""
        maf = self.df
        accession = maf.get_value(gene, organism)
        return accession

    def get_accessions(self, go_list=None):
        """Can take a gene/organism list as an argument:
                go_list = [[gene.1, org.1], ... , [gene.n, org.n]]
        Or it takes an empty gene/organisms list, which returns the 
        entire list of accession numbers."""
        if go_list is None:
            go_list = self.go_list
        maf = self.df
        genes = []
        organisms = []
        for gene, organism in go_list:
            genes.append(gene)
            organisms.append(organism)
        accessions = maf.lookup(genes, organisms).tolist()
        return accessions

    def get_accession_alignment(self, gene):
        """Takes a single gene and returns a list of accession numbers
        for the different orthologs."""
        maf = self.df
        accession_alignment = maf.T[gene].tolist()[1:]
        return accession_alignment

    def get_tier_frame(self, tiers=None):
        """Takes a list of tiers or nothing.
        Returns a dictionary as:
            keys:  Tier list
            values:  Data-Frame for each tier"""
        maf = self.df
        tier_frame_dict = {}
        if tiers is None:
            tiers = maf.groupby('Tier').groups.keys()
        for tier in tiers:
            tier = str(tier)
            tier_frame_dict[tier] = maf.groupby('Tier').get_group(tier)
        return tier_frame_dict

    def get_acc_dict(self):
        """This function takes a list of accession numbers and returns a dictionary
        which contains the corresponding genes/organisms."""
        gene_list = self.gene_list
        org_list = self.org_list
        maf = self.df
        go1 = []
        go = {}
        for gene in gene_list:
            for org in org_list:
                query_acc = self.get_accession(gene, org)
                go1.append(gene)
                go1.append(org)
                go[query_acc] = go1
                go1 = []
        return go

    def get_master_lists(self, df):
        """This function populates the organism and gene lists with a data frame.
        This function also populates several dictionaries.
        The dictionaries contain separate keys for Missing genes."""

        maf = df
        self.gene_list = maf.index.tolist()
        self.gene_count = len(self.gene_list)

        self.org_list = maf.axes[1].tolist()[1:]
        self.org_count = len(self.org_list)

        self.tier_list = maf['Tier'].tolist()
        self.tier_dict = maf['Tier'].to_dict()
        self.tier_frame_dict = self.get_tier_frame()

        self.acc_dict = self.get_acc_dict()
        self.acc_list = list(self.acc_dict.keys())

