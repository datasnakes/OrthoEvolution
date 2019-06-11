"""Comparative Genetics"""
# Standard Library
import os
import shutil
import time
import copy
import random
import pkg_resources
from pathlib import Path
# OrthoEvol
from OrthoEvol.Manager.config import data
from OrthoEvol.Manager.management import ProjectManagement
from OrthoEvol.utilities import FullUtilities
from OrthoEvol.Tools.logit import LogIt
# Other
import pandas as pd
from ete3 import NCBITaxa
# from pandas import ExcelWriter
# NCBITaxa().update_taxonomy_database()

# TODO: Create function for archiving and multiple runs (this can go
# into the Management class)


class BaseComparativeGenetics(object):
    """Base class in the Blast module."""

    __acc_filename = ''
    __paml_filename = ''
    __acc_path = ''
    __data = ''

    # Initialize Logging
    blastn_log = LogIt().default(logname="blastn", logfile=None)

    # TODO:  CREAT PRE-BLAST and POST-BLAST functions
    def __init__(self, project=None, project_path=os.getcwd(), acc_file=None,
                 taxon_file=None, pre_blast=False, post_blast=True, hgnc=False,
                 proj_mana=None, **kwargs):
        """This is the base class for the Blast module.

        It parses an accession file in order to provide easy handling for data.

        The .csv accession file contains the following header info:
            * "Tier" - User defined.
            * "Gene" - HUGO Gene Nomenclature Committee(HGNC) symbol for the genes of interest.
            * Query Organism - A well annotated query organism.
            * Other organisms - The other headers are Genus_species of other taxa.

        The organisms are taken from:
        ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/multiprocessing/
        The genes are taken from:
        http://www.guidetopharmacology.org/targets.jsp.
        The API gives the user access to their data in a higher level for
        downstream processing or for basic observation of the data.

        :param project:  The name of the project.
        :param project_path:  The location of the project, which is generally
                            defined by the ProjectManagement configuration.
        :param acc_file:  The name of the accession file.
        :param taxon_file:  A file that contains an ordered list of taxonomy ids.
        :param pre_blast:  A flag that gives the user access to an API that
                        contains extra information about their genes using the
                        mygene package.
        :param post_blast:  A flag that is used to handle a BLAST result file,
                            which returns information about misssing
                            data, duplicates, etc.
        :param hgnc:  A flag used as a placeholder for future work with HGNC files.
        :param proj_mana:  This parameter is used to compose (vs inherit) the
                           ProjectManagement class with the ComparativeGenetics class.
                           This parameter allows the various blast classes to function with or
                           without the Manager module.
        :param kwargs:  The kwargs here are generally used for standalone blasting or for development.
        :returns:  A pandas data-frame, pivot-table, and associated lists and dictionaries.
        """

        # Private Variables
        self.__pre_blast = pre_blast
        self.__post_blast = post_blast
        self.__taxon_filename = taxon_file
        self.acc_file = acc_file
        self.acc_csv_filename = acc_file

        # Initialize variables
        self.project_path = project_path
        self.project = project
        self.species = 'Homo_sapiens'
        self.taxon_file = taxon_file
        self.go_list = None
        self.proj_mana = proj_mana

        # Initialize Logging
        self.get_time = time.time
        self.sep = 50*'*'

        # Initialize Utilities
        self.blast_utils = FullUtilities()
        if self.project_path and self.project:
            self.project_path = Path(project_path) / Path(project)
        elif self.project and not self.project_path:
            self.project_path = self.project
        # If user does not want to use a project name, create one anyway.
        elif not self.project and not self.project_path:
            four_ints = random.sample(range(1, 9), 4)  # 4 random integers
            four_ints_str = ''.join(str(e) for e in four_ints)
            self.project = "orthoevol" + four_ints_str  # New project name
            self.project_path = os.getcwd()
        # TODO: Add ability to use an existing project or project path.
        elif self.project_path and not self.project:
            raise NotImplementedError

        self.blastn_log.debug('Project name: %s' % self.project)
        self.blastn_log.debug('Project path: %s' % self.project_path)

        # Configuration of class attributes.
        add_self = self.blast_utils.attribute_config(cls=self,
                                                     composer=proj_mana,
                                                     checker=ProjectManagement,
                                                     project=project,
                                                     project_path=project_path)
        for variable, attribute in add_self.__dict__.items():
            setattr(self, variable, attribute)

        # Handle the taxon_id file and blast query
        if self.taxon_file is not None:
            # File init
            self.taxon_path = self.project_index / Path(self.__taxon_filename)
        # Handle the master accession file (could be before or after blast)
        if self.copy_from_package:
            shutil.copy(pkg_resources.resource_filename(data.__name__, self.acc_file),
                        str(self.project_index))
        else:
            shutil.copy(self.acc_file, str(self.project_index))
        self.acc_filename = self.acc_file
        if self.acc_file is not None:
            # File init
            self.acc_sqlite_filename = Path(acc_file).stem + '.sqlite'
            self.acc_sqlite_tablename = Path(acc_file).stem.replace('.', '_')
            self.acc_csv_path = self.project_index / Path(self.acc_csv_filename)
            self.acc_sqlite_path = self.project_index / Path(self.acc_sqlite_filename)

            # Handles for organism lists
            self.org_list = []
            self.ncbi_orgs = []
            self.org_count = 0
            self.taxon_ids = []
            self.taxon_orgs = []
            self.taxon_dict = {}

            # Handles for gene lists
            self.gene_list = []
            self.gene_count = 0

            # Handles for tier lists
            self.tier_list = []
            self.tier_dict = {}
            self.tier_frame_dict = {}

            # Handles for accession lists
            self.acc_dict = {}
            self.acc_list = []

            # Handles for blast queries
            self.blast_human = []
            self.blast_rhesus = []

            # Handles for different dataframe initializations
            self.raw_acc_data = self.blast_utils.accession_sqlite2pandas(self.acc_sqlite_tablename, self.acc_sqlite_filename,
                                                                                  path=self.project_index, acc_file=self.acc_csv_filename)
            # Master accession file for the blast
            self.building_filename = str(self.acc_file[:-4] + 'building.csv')

            # Pre-Blast objects
            self.mygene_df = pd.DataFrame()  # MyGene
            self.mygene_filename = "%s_mygene.csv" % self.project  # MyGene
            self.mygene_path = self.data / Path(self.mygene_filename)  # MyGene
            self.header = self.raw_acc_data.axes[1].tolist()

            # Blast accession numbers
            self.building = copy.deepcopy(self.raw_acc_data)
            del self.building['Tier']
            del self.building[self.species]
            # Object for good user output
            self.building = self.building.set_index('Gene')
            self.building_file_path = self.data / Path(self.building_filename)

            # Blast time points
            # Master time file for the blast
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')  # Master time file for the blast
            self.building_time = copy.deepcopy(self.raw_acc_data)
            del self.building_time['Tier']
            del self.building_time[self.species]
            self.building_time = self.building_time.set_index('Gene')
            self.building_time_file_path = self.data / Path(self.building_time_filename)

            # Handles for accession file analysis # #
            if self.__post_blast:
                # Missing
                self.missing_dict = {}
                self.missing_genes = {}
                self.missing_organsims = {}
                self.missing_gene_count = 0
                self.missing_organsims_count = 0

                # Duplicates
                self.duplicated_dict = {}
                self.duplicated_accessions = {}
                self.dup_acc_count = {}
                self.duplicated_genes = {}
                self.dup_gene_count = {}
                self.duplicated_organisms = {}
                self.dup_org_count = {}
                self.duplicated_random = {}
                self.duplicated_other = {}
                self.time_dict = {}

            # Format the main data frame #### #
            self.__data = self.raw_acc_data.set_index('Gene')
            self.df = self.__data
            # Format the main pivot table #### #
            self.pt = pd.pivot_table(
                copy.deepcopy(self.raw_acc_data),
                index=['Tier', 'Gene'],
                aggfunc='first')
            array = self.pt.axes[1].tolist()  # Organism list
            self.pt.columns = pd.Index(array, name='Organism')

            # Handles for full dictionaries #### #
            self.org_dict = self.df.ix[0:, self.species:].to_dict()
            self.gene_dict = self.df.T.to_dict()
            self.get_master_lists(self.__data)  # populates our lists
        else:
            self.building_filename = str(self.project + '_building.csv')
            self.building_time_filename = str(self.project + '_building_time.csv')


# //TODO-ROB Add HGNC python module
    @staticmethod
    def get_file_list(file):
        """Turn csv column to list.

        :param file: Name of csv file.
        """
        file_data = pd.read_csv(file, header=None)
        file_list = list(file_data[0])
        return file_list

    def get_master_lists(self, df, csv_file=None):
        """Populate the organism and gene lists with a data frame.

        It will also populate pre-blast attributes (mygene) and post-blast
        attributes (missing and duplicates) under the proper conditions.

        :param df: The preferred way of utilizing the function is with a data-frame.
        :param csv_file: If a csv_file is given, then a data-frame will be
                         created by reinitializing the object.
                         (Default value = None)
        :returns:  An API can be utilized to access a gene list, organism list,
                   taxon-id list, tier list/dict/data-frame, accession
                   list/data-frame, blast query list, mygene information, and
                   missing/duplicate information.
        """

        # Usually only a user would manually add a csv file for their own
        # purposes.
        self.blastn_log.info("Getting the master lists.")
        if csv_file is not None:
            self.__init__(project=self.project, acc_file=csv_file)
            df = self.df
        maf = df
        self.gene_list = maf.index.tolist()
        self.gene_count = len(self.gene_list)

        self.org_list = maf.axes[1].tolist()[1:]
        self.org_count = len(self.org_list)
        self.ncbi_orgs = list(org.replace('_', ' ') for org in self.org_list)

        if self.__taxon_filename is not None:
            # Load taxon ids from a file
            self.taxon_ids = self.get_file_list(self.taxon_path)
        else:
            # Load taxon ids from a local NCBI taxon database via ete3
            ncbi = NCBITaxa()
            taxon_dict = ncbi.get_name_translator(self.ncbi_orgs)
            self.taxon_ids = list(tid[0] for tid in taxon_dict.values())
            self.taxon_orgs = list(torg for torg in taxon_dict.keys())
            self.taxon_orgs = list(org.replace(' ', '_')
                                   for org in self.taxon_orgs)
            self.taxon_dict = dict(zip(self.taxon_orgs, self.taxon_ids))
            self.taxon_lineage = self.get_taxon_dict()

        self.tier_list = maf['Tier'].tolist()
        self.tier_dict = maf['Tier'].to_dict()
        self.tier_frame_dict = self.get_tier_frame()

        self.acc_dict = self.get_acc_dict()
        self.acc_list = list(self.acc_dict.keys())

        # Get blast query list
        if self.species == 'Homo_sapiens':
            self.blast_human = self.df.Homo_sapiens.tolist()
            self.blast_rhesus = self.df.Macaca_mulatta.tolist()

        # Pre-Blast gene analysis
        if self.__pre_blast is True:
            self.mygene_df = self.blast_utils.my_gene_info(
                acc_dataframe=copy.deepcopy(self.raw_acc_data))
            self.mygene_df.to_csv(self.mygene_path, index=False)

        # Post-Blast accession analysis
        if self.__post_blast:

            # Missing
            self.missing_dict = self.blast_utils.get_miss_acc(
                acc_dataframe=copy.deepcopy(self.raw_acc_data))
            self.missing_genes = self.missing_dict['genes']
            self.missing_gene_count = self.missing_genes['count']
            del self.missing_genes['count']
            self.missing_organsims = self.missing_dict['organisms']
            self.missing_organsims_count = self.missing_organsims['count']
            del self.missing_organsims['count']

            # Duplicates
            self.duplicated_dict = self.blast_utils.get_dup_acc(self.acc_dict, self.gene_list,
                                                                self.org_list)
            self.duplicated_accessions = self.duplicated_dict['accessions']
            self.duplicated_organisms = self.duplicated_dict['organisms']
            self.duplicated_genes = self.duplicated_dict['genes']
            self.duplicated_random = self.duplicated_dict['random']
            self.duplicated_other = self.duplicated_dict['other']

    def get_accession(self, gene, organism):
        """Access a single accession number.

        :param gene:  An input gene.
        :param organism:  An input organism.
        :return:  A single accession number of the target gene/organism.
        """

        maf = self.df
        accession = maf.at[gene, organism]
        if isinstance(accession, float):
            accession = 'missing'
        return accession

    def get_orthologous_gene_sets(self, go_list=None):
        """Access a list of accession numbers.

        :param go_list:  A nested list of gene/organism lists
                         (go_list = [[gene.1, org.1], ... , [gene.n, org.n]]).
                         (Default value = None)
        :return:  An ordered list of accession numbers (or "missing") that
                  correspond to the go_list index.
        """

        if self.go_list is None:
            accessions = self.acc_list
        else:
            accessions = []
            for gene, organism in self.go_list:
                accession = self.get_accession(gene, organism)
                accessions.append(accession)
        return accessions

    def get_orthologous_accessions(self, gene):
        """Take a single gene & return a list of accession numbers for the different orthologs.

        :param gene:  An input gene from the accession file.
        :return:  A list of accession numbers that correspond to the orthologs
                  of the target gene.
        """

        maf = self.df
        accession_alignment = maf.T[gene].tolist()[1:]
        return accession_alignment

    def get_tier_frame(self, tiers=None):
        """Organize a dictionary by tier.

        Each tier (key) has a value, which is a data-frame of genes
        associated with that tier.

        :param tiers:  A list of tiers in the accession file.
                       (Default value = None)
        :return:  A nested dictionary for accessing information by tier.
        """

        maf = self.df
        tier_frame_dict = {}
        if tiers is None:
            tiers = maf.groupby('Tier').groups.keys()
        for tier in tiers:
            tier_frame_dict[str(tier)] = maf.groupby('Tier').get_group(tier)
        return tier_frame_dict

    def get_taxon_dict(self):
        """Get the taxonomy information about each organism using ETE3.

        :return:  Returns several dictionaries.  One is a basic organism (key)
        to taxonomy id (value) dictionary, and the other is a lineage
        dictionary with the an organism key and a lineage dictionary as the
        value.  The lineage dictionary keys for each organism are
        ["class", "family", "genus", "kingdowm", "order", "phylum", "species",
        "superkingdom"].
        """

        ncbi = NCBITaxa()
        taxa_dict = {}
        for organism in self.org_list:
            taxa_dict[organism] = {}
            taxid = self.taxon_dict[organism]
            lineage = ncbi.get_lineage(taxid)
            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            for id in lineage:
                if ranks[id] == 'no rank':
                    continue
                if ranks[id] not in ['superkingdom', 'kingdom', 'phylum',
                                     'class', 'order', 'family', 'genus', 'species']:
                    continue
                taxa_dict[organism][ranks[id]] = names[id]
        return taxa_dict

    def get_acc_dict(self):
        """Input a list of accession numbers and return a dictionary with corresponding genes/organisms.

        :return: An accession dictionary who's values are nest gene/organism lists.
        """
        # TODO-ROB set up function to accept a parameter for unique values or
        # potential duplicates

        gene_list = self.gene_list
        org_list = self.org_list
        go = {}
        for gene in gene_list:
            for org in org_list:
                query_acc = self.get_accession(gene, org)
                if query_acc not in go:
                    go[query_acc] = []
                # TODO-ROB: Rework the missing functino using this.. maybe??
                elif query_acc == 'missing':
                    continue
                self.go_list = [gene, org]
                # Append so that duplicates can be identified
                go[query_acc].append(self.go_list)
        return go


class ComparativeGenetics(BaseComparativeGenetics):
    """Main Comparative Genetics class."""

    def __init__(self, project, template=None, taxon_file=None, post_blast=False, save_data=True, **kwargs):
        """Inherits BaseComparativeGenetics to build a file layer to the Blast workflow.

        This class handles all of the files before and after the Blast occurs.
        It also uses a building file to start where a previous blast left off.

        :param project:  The name of the project.
        :param template:  A template accession file in the desired format.
                          See the Blast README for an example.
        :param taxon_file:  A list of taxon ids in a text file.
        :param post_blast:  A flag that triggers the post blast analysis.
        :param save_data:  A flag that indicates whether the data should be
                           saved in an excel file or not.
        :param kwargs:  Mostly used for BaseComparativeGenetics
        :returns:  An API for accessing the various files used before, during,
                   and after blasting."""

        super().__init__(project=project, taxon_file=taxon_file,
                         post_blast=post_blast, hgnc=False, **kwargs)

        self.postblastlog = LogIt().default(logname="post blast", logfile=None)

        self.acc_file = template
        # Private variables
        self.__home = os.getcwd()
        if self.taxon_file is not None:
            self.__taxon_filename = self.taxon_file
            self.taxon_path = self.project_index / Path(self.taxon_file)
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
        """Build an accession file after a Blastn run.

        It finds whether or not the Blast has been interrupted or not, so that
        the Blast can pick up where it left off.

        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :param accession:  The accession of interest.
        :return:
        """

        # TODO-ROB:  Create this in the log file
        if pd.isnull(self.building.at[gene, organism]) is False:
            existing = self.building.at[gene, organism]
            if existing == accession:
                self.blastn_log.warning(self.sep)
                self.blastn_log.warning("Blastn has run on this gene.")
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
                self.blastn_log.warning("The existing accession will be overwritten.")

        self.building.at[gene, organism] = accession
        temp = self.building.reset_index()
        temp.insert(0, 'Tier', pd.Series(self.df['Tier'].tolist()))
        # TODO: make the query organism insert implicit
        temp.insert(2, self.species, self.df[self.species])
        temp.set_index('Tier')
        if self.save_data is True:
            temp.to_csv(str(self.building_file_path))

    def add_blast_time(self, gene, organism, start, end):
        """Build a file that stores the amount of time for each gene to blast.

        This method is similar to the add_accession() method.

        :param gene:  The gene of interest.
        :param organism:  The organism of interest.
        :param start:  Starting time.
        :param end:  Ending time.
        """
        # TODO-ROB Add a method that adds the time to the post-blast analysis API.
        # This will help us see if there is a correlation between gene, organism,
        # or accession with the length of time.

        elapsed_time = end - start
        # Edit the data frame
        self.building_time.at[gene, organism] = elapsed_time
        temp = self.building_time.reset_index()
        temp.insert(0, 'Tier', pd.Series(self.df['Tier'].tolist()))
        temp.insert(2, self.species, self.df[self.species])
        temp.set_index('Tier')
        if self.save_data is True:
            temp.to_csv(str(self.building_time_file_path))

    def post_blast_analysis(self, removed_genes=None):
        """Save the post blast data (duplicate/missing/removed) to an excel file.

        :param removed_genes:  Genes to exclude from the file.
                               (Default value = None)
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
            removed_worksheet = pd.DataFrame.from_dict(removed_genes_dict,
                                                       orient='index')
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
            msg = 'Random duplicates were added to your excel file.'
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
            gene_ms_count = pd.DataFrame.from_dict(gene_org_ms_count,
                                                   orient='index')
            gene_ms_count.to_excel(pb_file,
                                   sheet_name="Missing Organisms Count")
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
