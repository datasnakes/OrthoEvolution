import os
import shutil
from pathlib import Path

# NCBITaxa().update_taxonomy_database()
import pandas as pd
import pkg_resources
from ete3 import NCBITaxa

from Datasnakes.Manager import config
# from pandas import ExcelWriter
from Datasnakes.Manager.management import ProjectManagement
from Datasnakes.Orthologs.utils import attribute_config
from Datasnakes.Orthologs.Blast.utils import (my_gene_info, get_dup_acc,
                                              get_miss_acc)


# TODO-ROB Create function for archiving and multiple runs (this can go
# into the Management class)


class CompGenObjects(object):
    """ Comparative Genetics Analysis.

    Parses an accession file with the designated format in order to
    provide easy handling for data.  Creates python objects from the given
    data.

    Input:  An open .csv file object that contains a header of organisms.  The
    first column ranks the gene by tier, the second column is a HUGO Gene
    Nomenclature Committee(HGNC) symbol for the genes of interest.  The .csv
    has to be located in the same directory as this module unless a full path is
    specified.

    The organisms are taken from
    ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/multiprocessing/
    and the genes are taken from http://www.guidetopharmacology.org/targets.jsp.

    Output:  A pandas Data-Frame, Pivot-Table, and associated lists and dictionaries.
    """
    __acc_filename = ''
    __paml_filename = ''
    __acc_path = ''
    __data = ''

    # TODO-ROB:  CREAT PRE-BLAST and POST-BLAST functions
    def __init__(self, project=None, project_path=os.getcwd(), acc_file=None, taxon_file=None, pre_blast=False, post_blast=True, hgnc=False,
                 proj_mana=ProjectManagement, **kwargs):

        # Private Variables
        self.__pre_blast = pre_blast
        self.__post_blast = post_blast
        self.__taxon_filename = taxon_file
        self.acc_filename = acc_file

        self.project = project
        if project_path and project:
            self.project_path = Path(project_path) / Path(project)

        # Configuration of class attributes.
        add_self = attribute_config(self, composer=proj_mana, checker=ProjectManagement, project=project, project_path=project_path)
        for var, attr in add_self.__dict__.items():
            setattr(self, var, attr)

        # Handle the taxon_id file and blast query
        if taxon_file is not None:
            # File init
            self.taxon_path = self.project_index / Path(self.__taxon_filename)
        # Handle the master accession file (could be before or after blast)
        if kwargs['copy_from_package']:
            shutil.copy(pkg_resources.resource_filename(config.__name__, kwargs['MAF']), str(self.project_index))
            acc_file = kwargs['MAF']
            self.acc_filename = acc_file
        if acc_file is not None:
            # File init
            self.acc_path = self.project_index / Path(self.acc_filename)
            # Handles for organism lists #
            self.org_list = []
            self.ncbi_orgs = []
            self.org_count = 0
            self.taxon_ids = []
            self.taxon_orgs = []
            self.taxon_dict = {}
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
            # Handles for blast queries #
            self.blast_human = []
            self.blast_rhesus = []
            # Handles for different dataframe initializations#
            self.raw_acc_data = pd.read_csv(str(self.acc_path), dtype=str)
            self.building_filename = str(acc_file[:-4] + 'building.csv')  # Master accession file for the blast
            # #### Pre-Blast objects
            self.mygene_df = pd.DataFrame()  # MyGene
            self.mygene_filename = "%s_mygene.csv" % self.project  # MyGene
            self.mygene_path = self.data / Path(self.mygene_filename)  # MyGene
            self.header = self.raw_acc_data.axes[1].tolist()
            # #### Blast accession numbers
            self.building = pd.read_csv(str(self.acc_path), dtype=str)
            del self.building['Tier']
            del self.building['Homo_sapiens']
            self.building = self.building.set_index('Gene')  # Object for good user output
            self.building_file_path = self.data / Path(self.building_filename)
            # #### Blast time points
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')  # Master time file for the blast
            self.building_time = pd.read_csv(str(self.acc_path), dtype=str)
            del self.building_time['Tier']
            del self.building_time['Homo_sapiens']
            self.building_time = self.building_time.set_index('Gene')
            self.building_time_file_path = self.data / Path(self.building_time_filename)

            # # Handles for accession file analysis # #
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

            # #### Format the main data frame #### #
            self.__data = self.raw_acc_data.set_index('Gene')
            self.df = self.__data
            # #### Format the main pivot table #### #
            self.pt = pd.pivot_table(
                pd.read_csv(self.acc_path),
                index=['Tier','Gene'],
                aggfunc='first')
            array = self.pt.axes[1].tolist()  # Organism list
            self.pt.columns = pd.Index(array, name='Organism')

            # #### Handles for full dictionaries #### #
            self.org_dict = self.df.ix[0:, 'Homo_sapiens':].to_dict()
            self.gene_dict = self.df.T.to_dict()
            self.get_master_lists(self.__data)  # populates our lists
        else:
            self.building_filename = str(self.project + 'building.csv')
            self.building_time_filename = self.building_filename.replace('building.csv', 'building_time.csv')


# //TODO-ROB Add HGNC python module
    @staticmethod
    def get_file_list(file):
        data = pd.read_csv(file, header=None)
        file_list = list(data[0])
        return file_list

    def get_master_lists(self, df=None, csv_file=None):
        """This function populates the organism and gene lists with a data frame.

        This function also populates several dictionaries.
        The dictionaries contain separate keys for Missing genes.
        """

        # Usually a only a user would manually add a csv file for their own purposes.
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
        self.blast_human = self.df.Homo_sapiens.tolist()
        self.blast_rhesus = self.df.Macaca_mulatta.tolist()

        # Pre-Blast gene analysis
        if self.__pre_blast is True:
            self.mygene_df = my_gene_info(self.acc_path)
            self.mygene_df.to_csv(self.mygene_path, index=False)

        # Post-Blast accession analysis
        if self.__post_blast:
            # Missing
            self.missing_dict = get_miss_acc(self.acc_path)
            self.missing_genes = self.missing_dict['genes']
            self.missing_gene_count = self.missing_genes['count']
            del self.missing_genes['count']
            self.missing_organsims = self.missing_dict['organisms']
            self.missing_organsims_count = self.missing_organsims['count']
            del self.missing_organsims['count']
            # Duplicates
            self.duplicated_dict = get_dup_acc(self.acc_dict, self.gene_list, self.org_list)
            self.duplicated_accessions = self.duplicated_dict['accessions']
            self.duplicated_organisms = self.duplicated_dict['organisms']
            self.duplicated_genes = self.duplicated_dict['genes']
            self.duplicated_random = self.duplicated_dict['random']
            self.duplicated_other = self.duplicated_dict['other']

    def get_accession(self, gene, organism):
        """Get the accession of a gene.
        Args:
            gene = input a gene
            organism = input an organism
        Returns:
            accession = the accession will be returned.
        """
        maf = self.df
        accession = maf.get_value(gene, organism)
        if isinstance(accession, float):
            accession = 'missing'
        return accession

    def get_orthologous_gene_sets(self, go_list=None):
        """ Get a list of acessions.
        Args:
            Can take a gene/organism list as an argument:
            go_list = [[gene.1, org.1], ... , [gene.n, org.n]]
            It also takes an empty gene/organisms list.
        Returns:
            It returns an entire list of accession numbers.
        """
        if go_list is None:
            accessions = self.acc_list
        else:
            accessions = []
            for gene, organism in go_list:
                accession = self.get_accession(gene, organism)
                accessions.append(accession)
        return accessions

    def get_orthologous_accessions(self, gene):
        """Takes a single gene and returns a list of accession numbers
        for the different orthologs.
        """
        maf = self.df
        accession_alignment = maf.T[gene].tolist()[1:]
        return accession_alignment

    def get_tier_frame(self, tiers=None):
        """
        Args:
            Takes a list of tiers or nothing.
        Returns:
            It returns a dictionary.
            keys:  Tier list
            values:  Data-Frame for each tier.
        """
        maf = self.df
        tier_frame_dict = {}
        if tiers is None:
            tiers = maf.groupby('Tier').groups.keys()
        for tier in tiers:
            tier = str(tier)
            tier_frame_dict[tier] = maf.groupby('Tier').get_group(tier)
        return tier_frame_dict

    def get_taxon_dict(self, min=False):
        """Get the taxonomny information about each organism using ETE3.
        Returns:
            taxa_dict = dict of organisms and their taxonomy ids.
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
        # TODO-ROB set up function to accept a parameter for unique values or
        # potential duplicates
        """
        This function takes a list of accession numbers and returns a dictionary
        which contains the corresponding genes/organisms.

        :return: A dictionary with values that are nested lists.  Unique accessions
        will have a nested list with a length of 1.  Accessions that are duplicated
        will have a nested list with a length > 1.
        """
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
                go_list = [gene, org]
                # Append so that duplicates can be identified
                go[query_acc].append(go_list)
        return go
