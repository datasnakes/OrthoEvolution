import os
#import mygene
from ete3 import NCBITaxa
# NCBITaxa().update_taxonomy_database()
import pandas as pd
from pathlib import Path
# from pandas import ExcelWriter
from Datasnakes.Manager.utils.mana import ProjMana as ProjectManagement
import shutil
import pkg_resources
from Datasnakes.Manager import index
# TODO-ROB Create function for archiving and multiple runs (this can go
# into the Mana class)


class CompGenAnalysis(object):
    """ Comparative Genetics Analysis.

    Parses an accession file with the designated format in order to
    provide easy handling for data.

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
    def __init__(self, project, project_path=None, acc_file=None, taxon_file=None, paml_file=None, go_list=None, post_blast=True, hgnc=False,
                 proj_mana=ProjectManagement, **kwargs):
        # Private Variables
        self.__post_blast = post_blast
        self.__taxon_filename = taxon_file
        self.__paml_filename = paml_file
        self.acc_filename = acc_file
        self.project = project

        # Initiate the project management variable
        if not isinstance(proj_mana, ProjectManagement):
            if project_path:
                self.project_path = Path(project_path) / Path(self.project)
            else:
                self.project_path = Path(os.getcwd()) / Path(self.project)
            Path.mkdir(self.project_path, parents=True, exist_ok=True)
            print('project_path=%s' % self.project_path)
            self.removed_pm_config(kwargs)
        else:
            setattr(proj_mana, 'project', project)
            for key, value in proj_mana.__dict__.items():
                setattr(self, key, value)
            print('project_path=%s' % self.project_path)


        # Handle the taxon_id file and blast query
        if taxon_file is not None:
            # File init
            self.taxon_path = self.project_index / Path(self.__taxon_filename)
        # Handle the paml organism file
        # TODO-ROB Deprecate paml_file
        if paml_file is not None:
            # File init
            self.paml_path = self.project_index / Path(self.__paml_filename)
            self.paml_org_list = []
        # Handle the master accession file (could be before or after blast)
        if kwargs['copy_from_package']:
            shutil.copy(pkg_resources.resource_filename(index.__name__, kwargs['MAF']), str(self.project_index))
            acc_file = kwargs['MAF']
            self.acc_filename = acc_file
        if acc_file is not None:
            # File init
            self.acc_path = self.project_index / Path(self.acc_filename)
            self.raw_acc_data = pd.read_csv(str(self.acc_path), dtype=str)
            self.go_list = go_list
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
            # Handles for dataframe init #
            self.raw_acc_data = pd.read_csv(str(self.acc_path), dtype=str)

            self.building_filename = str(acc_file[:-4] + 'building.csv')
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')
            self.building = pd.read_csv(str(self.acc_path), dtype=str)
            del self.building['Tier']
            del self.building['Homo_sapiens']
            self.building = self.building.set_index('Gene')
            self.building_file_path = self.raw_data / \
                                      Path(self.building_filename)

            self.building_time = pd.read_csv(str(self.acc_path), dtype=str)
            del self.building_time['Tier']
            del self.building_time['Homo_sapiens']
            self.building_time = self.building_time.set_index('Gene')
            self.building_time_file_path = self.raw_data / \
                                           Path(self.building_time_filename)
            # self.mygene_df = pd.DataFrame()
            # self.mygene_filename = "%s_mygene.csv" % self.project
            # self.mygene_path = self.data / Path(self.mygene_filename)
            self.header = self.raw_acc_data.axes[1].tolist()

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
                pd.read_csv(
                    self.acc_path),
                index=[
                    'Tier',
                    'Gene'],
                aggfunc='first')
            array = self.pt.axes[1].tolist()  # Organism list
            self.pt.columns = pd.Index(array, name='Organism')

            # #### Handles for full dictionaries #### #
            self.org_dict = self.df.ix[0:, 'Homo_sapiens':].to_dict()
            self.gene_dict = self.df.T.to_dict()
            self.get_master_lists(self.__data)  # populates our lists
        else:
            self.building_filename = str(self.project + 'building.csv')
            self.building_time_filename = self.building_filename.replace(
                'building.csv', 'building_time.csv')

    def removed_pm_config(self, kwargs):
        self.project_index = self.project_path / Path('index')
        self.raw_data = self.project_path / Path('raw_data')
        self.data = self.project_path / Path('data')

        for key, value in kwargs.items():
            setattr(self, key, value)

        Path.mkdir(self.project_index, exist_ok=True)
        Path.mkdir(self.raw_data, exist_ok=True)
        Path.mkdir(self.data, exist_ok=True)


# TODO-ROB Add HGNC python module
    @staticmethod
    def get_file_list(file):
        data = pd.read_csv(file, header=None)
        file_list = list(data[0])
        return file_list

# ***********************************************PRE BLAST ANALYSIS TOOLS********************************************* #
# ***********************************************PRE BLAST ANALYSIS TOOLS********************************************* #
#     def my_gene_info(self):
#         # TODO-ROB
#         # TODO-SDH Add custom mygene options
#         # Initialize variables and import my-gene search command
#         urls = []
#         df = self.raw_acc_data
#
#         mg = mygene.MyGeneInfo()
#         # Create a my-gene query handle to get the data
#         human = list(x.upper() for x in self.blast_human)
#         mygene_query = mg.querymany(human, scopes='refseq',
#                                     fields='symbol,name,entrezgene,summary',
#                                     species='human', returnall=True, as_dataframe=True,
#                                     size=1, verbose=True)
#         # TODO-ROB:  Logging here
#         # Turn my-gene queries into a data frame and then reset the index
#         mygene_query['out'].reset_index(level=0, inplace=True)
#         mg_df = pd.DataFrame(mygene_query['out'])
#         mg_df.drop(mg_df.columns[[1, 2, 6]], axis=1, inplace=True)
#         # Rename the columns
#         mg_df.rename(columns={'entrezgene': 'Entrez ID', 'summary':
#                              'Gene Summary', 'query': 'RefSeqRNA Accession', 'name': 'Gene Name'},
#                      inplace=True)
#
#         # Create NCBI links using a for loop and the Entrez IDs
#         urls = [('<a href="{0}">{0}</a>'.format('https://www.ncbi.nlm.nih.gov/gene/' + str(entrez_id)))
#                 for entrez_id in mg_df['Entrez ID']]
#
#         ncbi = pd.DataFrame(urls, columns=['NCBI Link'], dtype=str)
#         # Merge, sort, and return the my-gene data frame
#
#         hot_data = pd.concat([pd.Series(df.Tier, dtype=str), df.Gene, mg_df, ncbi], axis=1)
#         hot_data.rename(columns={'Gene': 'Gene Symbol'}, inplace=True)
#         hot_data = hot_data.sort_values(['Tier'], ascending=True)
#
#         return hot_data

    def get_master_lists(self, df=None, csv_file=None):
        """This function populates the organism and gene lists with a data frame.

        This function also populates several dictionaries.
        The dictionaries contain separate keys for Missing genes.
        """
        if csv_file is not None:
            self.__init__(csv_file)
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
        if self.__paml_filename is not None:
            self.paml_org_list = self.get_file_list(self.paml_path)
        else:
            self.paml_org_list = self.paml_org_formatter()

        self.tier_list = maf['Tier'].tolist()
        self.tier_dict = maf['Tier'].to_dict()
        self.tier_frame_dict = self.get_tier_frame()

        self.acc_dict = self.get_acc_dict()
        self.acc_list = list(self.acc_dict.keys())

        # Get blast query list
        self.blast_human = self.df.Homo_sapiens.tolist()
        self.blast_rhesus = self.df.Macaca_mulatta.tolist()

        # Gene analysis
        # self.mygene_df = self.my_gene_info()
        # self.mygene_df.to_csv(self.mygene_path, index=False)

        # Accession file analysis
        if self.__post_blast:
            # Missing
            self.missing_dict = self.get_miss_acc()
            self.missing_genes = self.missing_dict['genes']
            self.missing_gene_count = self.missing_genes['count']
            del self.missing_genes['count']
            self.missing_organsims = self.missing_dict['organisms']
            self.missing_organsims_count = self.missing_organsims['count']
            del self.missing_organsims['count']
            # Duplicates
            self.duplicated_dict = self.get_dup_acc()
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

    def get_accessions(self, go_list=None):
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

    def get_accession_alignment(self, gene):
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

    # TODO-ROB turn into a utility
    # XXX PAML no longer needs a format different than `Homo_sapiens`
    def paml_org_formatter(self):
        org_list = []
        for organism in self.org_list:
            genus, sep, species = organism.partition('_')
            org = ''.join([genus[0], sep, species[0:28]])
            org_list.append(org)
        return org_list

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

    def get_dup_acc(self):
        """Get duplicate accessions.
        This function is used to analyze an accession file post-BLAST.
        It uses the accession dictionary as a base.

        :return: A master duplication dictionary used to initialize the
        duplicate class variables.
        """
        duplicated_dict = dict()
        duplicated_dict['accessions'] = {}
        duplicated_dict['genes'] = {}
        duplicated_dict['organisms'] = {}
        duplicated_dict['random'] = {}
        duplicated_dict['other'] = {}
        acc_dict = self.acc_dict
        for accession, go_list in acc_dict.items():
            # Finding duplicates by using the length of the accession
            # dictionary
            if len(go_list) > 1:
                # dict['accessions']['XM_000000'] = [[g, o], [g, o]]
                duplicated_dict['accessions'][accession] = go_list
                genes, orgs = zip(*go_list)
                genes = list(genes)
                orgs = list(orgs)
                # Process the duplicates by categorizing and storing in a
                # dictionary
                for go in go_list:
                    g = go[0]
                    o = go[1]
                    # Initialize the dictionaries if they haven't already been
                    if g not in duplicated_dict['genes']:
                        duplicated_dict['genes'][g] = {}
                    if o not in duplicated_dict['organisms']:
                        duplicated_dict['organisms'][o] = {}
                # Categorize the different types of duplication
                    # Duplicates that persist across an organisms
                    if orgs.count(o) == len(go_list):
                        duplicated_dict['organisms'][o][accession] = genes
                        del duplicated_dict['genes'][g]
                        break
                    # Duplication across an organisms, but also somewhere else
                    elif orgs.count(o) != 1:
                        alt_genes = list(
                            gene for gene, org in go_list if org == o)
                        duplicated_dict['organisms'][o][accession] = alt_genes

                    # Duplicates that persist across a gene
                    if genes.count(g) == len(go_list):
                        duplicated_dict['genes'][g][accession] = orgs
                        del duplicated_dict['organisms'][o]
                        break
                    # Duplication across a gene, but also somewhere else
                    elif genes.count(g) != 1:
                        alt_orgs = list(
                            org for gene, org in go_list if gene == g)
                        duplicated_dict['genes'][g][accession] = alt_orgs

                    # This is the "somewhere else" if the duplication is random or not categorized
                        # The duplication is random
                    if genes.count(g) == 1 and orgs.count(o) == 1:
                        del duplicated_dict['organisms'][o]
                        del duplicated_dict['genes'][g]
                        if accession not in duplicated_dict['random']:
                            duplicated_dict['random'][accession] = []
                        duplicated_dict['random'][accession].append(go)
                        # There is another category of duplication that I'm missing
                        # TODO-ROB:  If an other exists throw a warning in the
                        # logs
                    else:
                        del duplicated_dict['organisms'][o]
                        del duplicated_dict['genes'][g]
                        if accession not in duplicated_dict['other']:
                            duplicated_dict['other'][accession] = []
                        duplicated_dict['other'][accession].append(go)
            # Duplicate Organism count dictionary
            dup_org = pd.DataFrame.from_dict(self.duplicated_organisms)
            for org in self.org_list:
                try:
                    self.dup_org_count[org] = dup_org[org].count()
                except KeyError:
                    self.dup_org_count[org] = 0
            # Duplicate Gene count dictionary
            dup_gene = pd.DataFrame.from_dict(self.duplicated_genes)
            for gene in self.gene_list:
                try:
                    self.dup_gene_count[gene] = dup_gene[gene].count()
                except KeyError:
                    self.dup_gene_count[gene] = 0
            # Duplicate Accession count dictionary
            for accn, go in duplicated_dict['accessions'].items():
                self.dup_acc_count[accn] = go.__len__()
        return duplicated_dict

    def get_miss_acc(self, acc_file=None):
        """
        This function is used to analyze an accession file post BLAST.
        It generates several files and dictionaries regarding missing accession numbers.

        :param acc_file: An accession file (post BLAST).
        :return: A dictionary with data about the missing accession numbers by Gene and
        by Organism.
        """
        # TODO-ROB: Add Entrez ID validation;  Get info from xml files???
        missing_dict = dict()
        missing_dict['organisms'] = {}
        missing_dict['genes'] = {}

        # Initialize the Data Frames
        if acc_file is None:
            miss_gene_df = self.df.isnull()
            miss_org_df = self.df.T.isnull()
        else:
            self.__init__(acc_file=acc_file)
            miss_gene_df = self.df.isnull()
            miss_org_df = self.df.T.isnull()

        # Get missing Accessions by Organism
        organism_dict = miss_org_df.sum(axis=1).to_dict()
        total_miss = 0
        for organism, miss in organism_dict.items():
            if miss != 0:
                missing_dict['organisms'][organism] = {}
                # Missing Gene dict {'HTR1A': True}
                missing_genes = miss_gene_df.ix[:, organism].to_dict()
                # Do a list comprehension to get a list of genes
                missing_dict['organisms'][organism]['missing genes'] = list(key for key, value in missing_genes.items()
                                                                            if value)  # Value is True for miss accns

                # Number of missing accessions per organism
                missing_dict['organisms'][organism]['count'] = miss
                total_miss += miss
        missing_dict['organisms']['count'] = total_miss

        # Get missing Accessions by Gene
        gene_dict = miss_gene_df.sum(axis=1).to_dict()
        total_miss = 0
        for gene, miss in gene_dict.items():
            if miss != 0:
                missing_orgs = miss_gene_df.T.ix[:, gene].to_dict()
                missing_dict['genes'][gene] = {}
                # Do a list comprehension to get a list of organisms
                missing_dict['genes'][gene]['missing organisms'] = list(key for key, value in missing_orgs.items()
                                                                        if value  # Value is True for missing accessions
                                                                        if key != 'Tier')  # Don't include 'Tier'

                # Number of missing accessions per gene
                missing_dict['genes'][gene]['count'] = miss
                total_miss += miss
        missing_dict['genes']['count'] = total_miss

        return missing_dict

    def get_pseudogenes(self):
        """ UNDER DEVELOPMENT!!!
        This subclass will denote which genes are sudogenes.
        """
        print(self.__doc__)
