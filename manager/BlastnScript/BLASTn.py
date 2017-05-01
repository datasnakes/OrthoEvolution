#!/usr/local/apps/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
File Name: blastntest.py
Description: This script inputs a list of organisms & genes. It tests a script.

@authors: Robert A. Gilmore & Shaurita D. Hutchins
Date Created: March 29, 2017
Project Name: Orthologs Project

Edited and updated for use with local/standalone NCBI BLAST 2.6.0
"""
# Modules used
import os
import time   # Used to delay when dealing with NCBI server errors
import csv
import subprocess
import pandas as pd  # Used for dealing with CSV files
from Bio import SearchIO  # Used for parsing and sorting XML files.
from Bio.Blast.Applications import NcbiblastnCommandline  # Used for Local Blasting.
# Modules as custom
import logging as log
from datetime import datetime as d
# Proprietary Modules
from manager.lister import Lister
from manager.logit.logit import LogIt
from manager.BLASTingTemplate import BLASTingTemplate as BT
# TODO-ROB: Find packages for script timing and analysis


class BLAST(BT):
    def __init__(self, template=None):
        """Inherit from the BLASTing Template."""
        super().__init__(template=template)

        # Manage Directories
        self.__home = os.getcwd()
        self.__output = 'data/blastn/blast-xml-output/'  # Output directory
        self.processed = 'data/processed/'  # Processed data directory
        os.makedirs(self.__output, exist_ok=True)  # only make blast-xml-output dir since blastn dir exists
        # # Initialize Logging
        # self.__blastn_log = LogIt.blastn()
        # df = LogIt()
        # self.__date_format = df.date_format
        # self.get_time = time.time  # To get the time use 'get_time()'

        # ------------------------------------------------------------------------------
        self.blastn_log.info("These are the organisms: " + str(self.org_list))
        self.blastn_log.info("These are the genes: " + str(self.gene_list))
        self.blastn_log.info("These are the taxonomy ids: " + str(self.taxon_ids))
        # ---------------------------------------------------------------------

    @staticmethod
    def map_func(hit):
        """The map function for formatting hit id's.
        This will be used later in the script."""
        hit.id1 = hit.id.split('|')[3]
        hit.id2 = hit.id.split('|')[1]
        hit.id = hit.id[:-2]
        return hit

    def blast_file_config(self, file):
        """This function configures different files for new BLASTS.
        It also helps recognize whether or not a BLAST was terminated
        in the middle of the dataset.  This removes the last line of
        the accession file if it is incomplete."""
        os.chdir(self.__output)  # Change or move to the output directory
        output_dir_list = os.listdir()  # Make a list of files
        if file in output_dir_list:
            with open(file, 'r') as fi:
                f = csv.reader(fi)
                count = -1
                for row in f:
                    count += 1
                    ending = row
                gene = ending[1]
                taxid = self.taxon_ids[count]
                org = self.org_list[(len(ending) - 2)]

                ncbi = str("""result_handle1 = NcbiblastnCommandline(query="temp.fasta", db="refseq_rna", strand="plus", 
                evalue=0.001, out="%s_%s.xml", outfmt=5, gilist=%s + "gi", max_target_seqs=10, task="blastn")"""
                           % (gene, org, taxid))
                self.blastn_log.warning("An incomplete accession file was produced from the previous BLAST,"
                                        "which was terminated midway through the procedure.")

                self.blastn_log.info("The last row looks like: \n\t%s\n\t%s\n" % (self.header, ending))
                self.blastn_log.info("The BLAST ended on the following query: \n%s" % ncbi)
                if len(ending) < len(self.header):
                    self.blastn_log.info("Restarting the BLAST for the previous gene...")
                    count = count - 1
                continued_gene_list = list(x for i, x in enumerate(self.gene_list, 1) if i > count)
                self.blasting(continued_gene_list, self.org_list)
        else:
            self.blastn_log.info("A new BLAST started at %s" % self.get_time()())
            self.blasting(self.gene_list, self.org_list)

    def blasting(self, genes, organisms):
        os.chdir(self.__output)
        start_time = time.time()  # Variable used to check the processing time
        self.blastn_log.info("------------------------------------------------------------------")
        self.blastn_log.info("The script name is %s" % os.path.basename(__file__))
        self.blastn_log.info("The script began on %s" % str(d.now().strftime(self.__date_format)))
        self.blastn_log.info("------------------------------------------------------------------")

        # Factor in the new gene list and parse the human querys
        index = (len(self.gene_list)- len(genes))
        for query in self.blast_human[index:]:
            self.blastn_log.info('Human Accession: %s' % query)
            for gene in genes:
                try:
                    os.mkdir(gene)
                    self.blastn_log.info("Directory Created: %s" % gene)
                    os.chdir(gene)
                    self.blastn_log.info("Current Directory: " + str(os.getcwd()))
                    self.blastn_log.info("\n")

                except FileExistsError:
                    self.blastn_log.info("Directory already exists: %s" % gene)
                    os.chdir(gene)
                    self.blastn_log.info("Current Directory: " + str(os.getcwd()))

            self.blastn_log.info("Gene of Interest: %s" % gene)
            self.blastn_log.info("Gene Tier: %s" % self.tier_dict[gene])

