"""Comparative Genetics Files"""
import os
from pathlib import Path
import pandas as pd

from Datasnakes.Orthologs.Blast.comparative_genetics_objects import CompGenObjects
from Datasnakes.Tools.logit import LogIt


class CompGenFiles(CompGenObjects):
    def __init__(self, project, template=None, taxon_file=None, post_blast=False, save_data=True, **kwargs):
        """Inherit CompGenObjects to build a file layer to the Blast workflow.

        This class handles all of the files before and after the Blast occurs.
        It also uses a building file to start where a previous blast left off.

        :param project:  The name of the project.
        :param template:  A template accession file in the desired format.  See the Blast README for an example.
        :param taxon_file:  A list of taxon ids in a text file.
        :param post_blast:  A flag that triggers the post blast analysis.
        :param save_data:  A flag that indicates whether the data should be saved in an excel file or not.
        :param kwargs:  Mostly used for CompGenObjects
        :returns:  An API for accessing the various files used before, during, and after Blasting.
        """
        super().__init__(project=project, acc_file=template, taxon_file=taxon_file, post_blast=post_blast, hgnc=False, **kwargs)

        self.postblastlog = LogIt().default(logname="post blast", logfile=None)

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

    def add_accession(self, gene, organism, accession):
        """This method builds an accession file after a Blastn run.

        It finds whether or not the Blast has been interrupted or not, so that
        the Blast can pick up where it left off.

        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :param accession:  The accession of interest.
        :return:
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
        # TODO-ROB Add a method that adds the time to the post-blast analysis API.
        # This will help us see if there is a correlation between gene, organism,
        # or accession with the length of time.
        """Build a file that stores the amount of time for each gene to blast.

        This method is similar to the add_accession() method.

        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :param start:  Starting time.
        :param end:  Ending time.
        :return:
        """
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
        """Save the post blast data (duplicate/missing/removed) to an excel file.

        :param removed_genes:  Genes to exclude from the file.
        :return:
        """
        # TODO-ROB  Fix the output format of the excel file.  View a sample
        # output in /Orthologs/comp_gen
        pba = '_postblastanalysis'
        pba_file_path = str(self.data / Path(self.project + pba + '.xlsx'))
        pb_file = pd.ExcelWriter(pba_file_path)

        # Removed Genes
        if removed_genes is not None:
            removed_genes_dict = {'Removed Genes': removed_genes}
            removed_worksheet = pd.DataFrame.from_dict(removed_genes_dict, orient='index')
            removed_worksheet.to_excel(pb_file, sheet_name="Removed Genes")
            msg = "Removed genes were added to your excel file."
            self.postblastlog.info(msg)

        # Duplicated Accessions
        try:
            acc_ws = pd.DataFrame.from_dict(self.dup_acc_count, orient='index')
            acc_ws.columns = ['Count']
            acc_ws.to_excel(pb_file, sheet_name="Duplicate Count by Accession")
            msg = "Dupilicate accessions were added to your excel file."
            self.postblastlog.info(msg)
        except (ValueError, AttributeError):
            pass

        # Duplicate Genes
        try:
            dup_gene_ws = pd.DataFrame.from_dict(
                self.dup_gene_count, orient='index')
            dup_gene_ws.columns = ['Count']
            dup_gene_ws.to_excel(pb_file, sheet_name="Duplicate Count by Gene")

            gene_org_dup = {}
            for gene, _ in self.duplicated_genes.items():
                gene_org_dup[gene] = []
                for _, genes in self.duplicated_genes[gene].items():
                    gene_org_dup[gene].append(genes)
            dup_org_ws2 = pd.DataFrame.from_dict(gene_org_dup, orient='index')
            dup_org_ws2.T.to_excel(
                pb_file, sheet_name="Duplicate Org Groups by Gene")
            msg = 'Dupilicate genes were added to your excel file.'
            self.postblastlog.info(msg)
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
            msg = 'Dupilicate species were added to your excel file.'
            self.postblastlog.info(msg)
        except (ValueError, AttributeError):
            pass

        # Random Duplicates
        try:
            rand_ws = pd.DataFrame.from_dict(
                self.duplicated_random, orient='index')
            rand_ws.to_excel(pb_file, sheet_name="Random Duplicates")
            msg = 'Random Duplicates were added to your excel file.'
            self.postblastlog.info(msg)
        except (ValueError, AttributeError):
            pass

        # Other Duplicates
        try:
            other_ws = pd.DataFrame.from_dict(
                self.duplicated_other, orient='index')
            other_ws.to_excel(pb_file, sheet_name="Other Duplicates")
            msg = 'Other duplicates were added to your excel file.'
            self.postblastlog.info(msg)
        except (ValueError, AttributeError):
            pass

        # Missing genes sorted by Organism
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

        # Missing Organisms sorted by Gene
        gene_org_ms = {}
        gene_org_ms_count = {}
        try:
            for gene, ms_dict in self.missing_genes.items():
                for key, value in ms_dict.items():
                    if key == 'missing genes':
                        gene_org_ms[gene] = value
                    else:
                        gene_org_ms_count[gene] = value
            gene_ms_count = pd.DataFrame.from_dict(gene_org_ms_count, orient='index')
            gene_ms_count.to_excel(pb_file, sheet_name="Missing Organisms Count")
            gene_ms = pd.DataFrame.from_dict(gene_org_ms, orient='index')
            gene_ms.to_excel(pb_file, sheet_name="Missing Organisms by Gene")
            msg = 'Missing Organisms by gene were added to your excel file.'
            self.postblastlog.exception(msg)
        except (ValueError, AttributeError):
            pass
        try:
            pb_file.save()
        except IndexError:
            msg = "There are no duplicates or missing genes."
            self.postblastlog.exception(msg)
