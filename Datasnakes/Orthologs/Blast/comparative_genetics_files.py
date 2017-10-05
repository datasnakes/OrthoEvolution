import os
import time
from pathlib import Path

import pandas as pd
from Datasnakes.Orthologs.Blast.comparative_genetics_objects import CompGenObjects
from Datasnakes.Tools import LogIt


# import pkg_resources
# import shutil


class CompGenFiles(CompGenObjects):
    """Perform Blast Analysis after completing CompGenBLASTn."""
    def __init__(self, project, template=None, taxon_file=None, post_blast=False, save_data=True, **kwargs):
        """Inherited from the CompGenObjects class.

        If the BLAST was cut short, then a build_file is to be used.
        """
        super().__init__(project=project, acc_file=template, taxon_file=taxon_file, post_blast=post_blast, hgnc=False, **kwargs)
        # TODO-ROB: Inherit or add variable for logger class
        # TODO-ROB Add Management directories
        # Private variables
        self.__home = os.getcwd()
        if taxon_file is not None:
            self.__taxon_filename = taxon_file
            self.taxon_path = self.project_index / Path(taxon_file)
        self.__post_blast = post_blast
        self.save_data = save_data

        if template is not None:
            self.template_filename = template
            self.template_path = self.project_index / \
                Path(self.template_filename)
            self.building_filename = str(template[:-4] + 'building.csv')
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')
        else:
            self.building_filename = str(self.project + 'building.csv')
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')

        # Initialize Logging
        logit = LogIt()
        self.blastn_log = logit.default('blastn', 'blastn.log')
        self.__date_format = logit._date_format
        self.get_time = time.time  # To get the time use 'get_time()'

        # Create variable for log separator
        log_sep = 50*'*'
        self.sep = log_sep

    def add_accession(self, gene, organism, accession):
        """Take an accession and add in to the building dataframe & csv file.

        This returns a log.
        """
        # TODO-ROB:  Create this in the log file
        if pd.isnull(self.building.get_value(gene, organism)) is False:
            existing = self.building.get_value(gene, organism)
            if existing == accession:
                self.blastn_log.warning(self.sep)
                self.blastn_log.warning("BlastN has run on this gene.")
                self.blastn_log.warning("The ACCESSION(%s) for the %s %s gene already exists in our data set."
                                        % (accession, organism, gene))
                self.blastn_log.warning(self.sep)
            else:
                self.blastn_log.critical(self.sep)
                self.blastn_log.critical("BlastN has run on this gene.")
                self.blastn_log.critical("The queried ACCESSION (%s) does not match the existing ACCESSION (%s)"
                                         % (accession, existing))

                self.blastn_log.critical("Queried Accession Gene: %s" % gene)
                self.blastn_log.critical(
                    "Queried Accession Organism: %s" %
                    organism)
                self.blastn_log.critical(
                    "Existing Accession Gene: %s" %
                    self.acc_dict[existing][0][0])
                self.blastn_log.critical(
                    "Existing Accession Organism: %s" %
                    self.acc_dict[existing][0][1])
                self.blastn_log.critical(self.sep)
                raise ValueError("The queried ACCESSION (%s) does not match the existing ACCESSION (%s).  Please see"
                                 "the log file." % (accession, existing))

        self.building.set_value(gene, organism, accession)
        temp = self.building.reset_index()
        temp.insert(0, 'Tier', pd.Series(self.df['Tier'].tolist()))
        # TODO-ROB make the query organism insert implicit
        temp.insert(2, 'Homo_sapiens', self.df['Homo_sapiens'])
        temp.set_index('Tier')
        if self.save_data is True:
            temp.to_csv(str(self.building_file_path))

    def add_blast_time(self, gene, organism, start, end):
        """Retrieve the start/end times and add in to the building_time dataframe & csv file."""
        elapsed_time = end - start
        # Edit the data frame
        self.building_time.set_value(gene, organism, elapsed_time)
        temp = self.building_time.reset_index()
        temp.insert(0, 'Tier', pd.Series(self.df['Tier'].tolist()))
        temp.insert(2, 'Homo_sapiens', self.df['Homo_sapiens'])
        temp.set_index('Tier')
        if self.save_data is True:
            temp.to_csv(str(self.building_time_file_path))

    def post_blast_analysis(self, removed_genes=None):
        """Analyze the blast results.

        Generate information about any duplicated or missing accessions by gene
        and by organism.
        """
        # TODO-ROB  Fix the output format of the excel file.  View a sample
        # output in /Orthologs/comp_gen
        pba_file_path = str(self.data / Path(self.project + '_pba.xlsx'))
        pb_file = pd.ExcelWriter(pba_file_path)

        # Removed Genes
        if removed_genes is not None:
            rm_ws = pd.DataFrame(removed_genes)
            rm_ws.to_excel(pb_file, sheet_name="Removed Genes")

        # Duplicated Accessions
        try:
            acc_ws = pd.DataFrame.from_dict(self.dup_acc_count, orient='index')
            acc_ws.columns = ['Count']
            acc_ws.to_excel(pb_file, sheet_name="Duplicate Count by Accession")
        except (ValueError, AttributeError):
            pass

        # Duplicate Genes
        try:
            dup_gene_ws = pd.DataFrame.from_dict(
                self.dup_gene_count, orient='index')
            dup_gene_ws.columns = ['Count']
            dup_gene_ws.to_excel(pb_file, sheet_name="Duplicate Count by Gene")

            gene_org_dup = {}
            for gene, dup_dict in self.duplicated_genes.items():
                gene_org_dup[gene] = []
                for acc, genes in self.duplicated_genes[gene].items():
                    gene_org_dup[gene].append(genes)
            dup_org_ws2 = pd.DataFrame.from_dict(gene_org_dup, orient='index')
            dup_org_ws2.T.to_excel(
                pb_file, sheet_name="Duplicate Org Groups by Gene")
        except (ValueError, AttributeError):
            pass

        # Species Duplicates
        try:
            dup_org_ws1 = pd.DataFrame.from_dict(
                self.dup_org_count, orient='index')
            dup_org_ws1.columns = ['Count']
            dup_org_ws1.to_excel(pb_file, sheet_name="Duplicate Count by Org")

            org_gene_dup = {}
            for gene, dup_dict in self.duplicated_organisms.items():
                org_gene_dup[gene] = []
                for acc, genes in self.duplicated_organisms[gene].items():
                    org_gene_dup[gene].append(genes)
            dup_org_ws2 = pd.DataFrame.from_dict(org_gene_dup, orient='index')
            dup_org_ws2.T.to_excel(
                pb_file, sheet_name="Duplicate Gene Groups by Org")
        except (ValueError, AttributeError):
            pass

        # Random Duplicates
        try:
            rand_ws = pd.DataFrame.from_dict(
                self.duplicated_random, orient='index')
            rand_ws.to_excel(pb_file, sheet_name="Random Duplicates")
        except (ValueError, AttributeError):
            pass

        # Other Duplicates
        try:
            other_ws = pd.DataFrame.from_dict(
                self.duplicated_other, orient='index')
            other_ws.to_excel(pb_file, sheet_name="Other Duplicates")
        except (ValueError, AttributeError):
            pass

        # Missing by Organism
        org_gene_ms = {}
        org_gene_ms_count = {}
        try:
            for org, ms_dict in self.missing_organsims.items():
                for key, value in ms_dict.items():
                    if key == 'missing genes':
                        org_gene_ms[org] = value
                    else:
                        org_gene_ms_count[org] = value
            org_ms_count = pd.DataFrame.from_dict(
                org_gene_ms_count, orient='index')
            org_ms_count.to_excel(pb_file, sheet_name="Missing Genes Count")
            org_ms = pd.DataFrame.from_dict(org_gene_ms, orient='index')
            org_ms.to_excel(pb_file, sheet_name="Missing Genes by Org")
        except (ValueError, AttributeError):
            pass

        # Missing by Gene
        gene_org_ms = {}
        gene_org_ms_count = {}
        try:
            for gene, ms_dict in self.missing_genes.items():
                for key, value in ms_dict.items():
                    if key == 'missing genes':
                        gene_org_ms[gene] = value
                    else:
                        gene_org_ms_count[gene] = value
            gene_ms_count = pd.DataFrame.from_dict(
                gene_org_ms_count, orient='index')
            gene_ms_count.to_excel(
                pb_file, sheet_name="Missing Organisms Count")
            gene_ms = pd.DataFrame.from_dict(gene_org_ms, orient='index')
            gene_ms.to_excel(pb_file, sheet_name="Missing Organisms by Genes")
        except (ValueError, AttributeError):
            pass
        try:
            pb_file.save()
        except IndexError:
            print("There are no duplicates or missing genes.")
