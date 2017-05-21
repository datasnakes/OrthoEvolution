import time
import os
from pathlib import Path
import pandas as pd
from Orthologs.CompGenetics.comp_gen import CompGenAnalysis as CGA
from Manager.logit.logit import LogIt


class BLASTAnalysis(CGA):
    def __init__(self, repo, user, project, research, research_type,
                 template=None, taxon_file=None, post_blast=False, save_data=True, **kwargs):
        """Inherit from the CompGenAnalysis class.  If the BLAST was cut short,
        then a build_file is to be used."""
        super().__init__(repo=repo, user=user, project=project, research=research, research_type=research_type,
                         acc_file=template, taxon_file=taxon_file, post_blast=post_blast, hgnc=False, **kwargs)
        # TODO-ROB: Inherit or add variable for logger class
        # TODO-ROB Add Mana directories
        # Private variables
        self.__home = os.getcwd()
        if taxon_file is not None:
            self.__taxon_filename = taxon_file
            self.taxon_path = self.project_index / Path(taxon_file)
        self.__post_blast = post_blast
        self.save_data = save_data

        if template is not None:
            self.template_filename = template
            self.template_path = self.project_index / Path(self.template_filename)
            self.building_filename = str(template[:-4] + 'building.csv')
            self.building_time_filename = self.building_filename.replace('building.csv', 'building_time.csv')
        else:
            self.building_filename = str(project + 'building.csv')
            self.building_time_filename = self.building_filename.replace('building.csv', 'building_time.csv')

        # Initialize a data frame and file to add accession numbers to
        # Initialize a data frame and file to add blast times to
        # self.building = self.raw_acc_data
        # del self.building['Tier']
        # del self.building['Homo_sapiens']
        # self.building = self.building.set_index('Gene')
        # self.building_time = self.building
        # self.building_file_path = self.raw_data / Path(self.building_filename)
        # self.building_time_file_path = self.raw_data / Path(self.building_time_filename)

        # Initialize Logging
        df = LogIt('blast_test.log', 'blastn')
        self.blastn_log = df.basic
        #self.postblast_log = df.basic
        #self.config_log = df.basic(self.user_log / Path('BLAST.log'))
        self.__date_format = df.date_format
        self.get_time = time.time  # To get the time use 'get_time()'
        # Logging variables
        # self.__date_format = '%a %b %d at %I:%M:%S %p %Y'  # Used to add as a date
        # self.__archive_format = '%m-%d-%Y@%I:%M:%S-%p'  # Used to append to archives
        # self.__log_format = '%(name)s - [%(levelname)-2s]: %(message)s'
        # log.basicConfig(level=log.DEBUG,
        #                 format=self.__log_format,
        #                 filename="logs/accessions2blastxml_%s.log" % str(d.now().strftime(self.__archive_format)))
        # self.blast_log = log.getLogger('Blastn')

    def add_accession(self, gene, organism, accession):
        """Takes an accession and adds in to the building dataframe,
        and also adds to the csv file.  This returns a log."""
        # TODO-ROB:  Create this in the log file
        if pd.isnull(self.building.get_value(gene, organism)) is False:
            existing = self.building.get_value(gene, organism)
            if existing == accession:
                self.blastn_log.warning("***********************************************************************")
                self.blastn_log.warning("This gene has already been BLASTed...")
                self.blastn_log.warning("The ACCESSION(%s) for the %s %s gene already exists in our data set."
                                        % (accession, organism, gene))
                self.blastn_log.warning("***********************************************************************")
            else:
                self.blastn_log.critical("***********************************************************************")
                self.blastn_log.critical("The gene has already been BLASTed...")
                self.blastn_log.critical("The queried ACCESSION (%s) does not match the existing ACCESSION (%s)"
                                         % (accession, existing))

                self.blastn_log.critical("Queried Accession Gene: %s" % gene)
                self.blastn_log.critical("Queried Accession Organism: %s" % organism)
                self.blastn_log.critical("Existing Accession Gene: %s" % self.acc_dict[existing][0][0])
                self.blastn_log.critical("Existing Accession Organism: %s" % self.acc_dict[existing][0][1])
                self.blastn_log.critical("***********************************************************************")
                raise ValueError("The queried ACCESSION (%s) does not match the existing ACCESSION (%s).  Please see"
                                 "the log file." % (accession, existing))

        self.building.set_value(gene, organism, accession)
        temp = self.building.reset_index()
        temp.insert(0, 'Tier', self.df['Tier'])
        if self.save_data is True:
            temp.to_csv(str(self.building_file_path))

    def add_blast_time(self, gene, organism, start, end):
        """Takes the start/end times and adds in to the building_time 
        dataframe, and also adds to the csv file."""
        elapsed_time = end - start
        # Edit the data frame
        self.building_time.set_value(gene, organism, elapsed_time)
        temp = self.building_time.reset_index()
        temp.insert(0, 'Tier', self.df['Tier'])
        if self.save_data is True:
            temp.to_csv(str(self.building_time_file_path))

    def post_blast_analysis(self, project_name, removed_genes=None):
        # TODO-ROB  Fix the output format of the excel file.  View a sample output in /Orthologs/comp_gen
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
        except ValueError:
            pass
        # Duplicate Genes
        try:
            dup_gene_ws = pd.DataFrame.from_dict(self.dup_gene_count, orient='index')
            dup_gene_ws.columns = ['Count']
            dup_gene_ws.to_excel(pb_file, sheet_name="Duplicate Count by Gene")

            gene_org_dup = {}
            for gene, dup_dict in self.duplicated_genes.items():
                gene_org_dup[gene] = []
                for acc, genes in self.duplicated_genes[gene].items():
                    gene_org_dup[gene].append(genes)
            dup_org_ws2 = pd.DataFrame.from_dict(gene_org_dup, orient='index')
            dup_org_ws2.T.to_excel(pb_file, sheet_name="Duplicate Org Groups by Gene")
        except ValueError:
            pass
        # Species Duplicates
        try:
            dup_org_ws1 = pd.DataFrame.from_dict(self.dup_org_count, orient='index')
            dup_org_ws1.columns = ['Count']
            dup_org_ws1.to_excel(pb_file, sheet_name="Duplicate Count by Org")

            org_gene_dup = {}
            for gene, dup_dict in self.duplicated_organisms.items():
                org_gene_dup[gene] = []
                for acc, genes in self.duplicated_organisms[gene].items():
                    org_gene_dup[gene].append(genes)
            dup_org_ws2 = pd.DataFrame.from_dict(org_gene_dup, orient='index')
            dup_org_ws2.T.to_excel(pb_file, sheet_name="Duplicate Gene Groups by Org")
        except ValueError:
            pass
        # Random Duplicates
        try:
            rand_ws = pd.DataFrame.from_dict(self.duplicated_random, orient='index')
            rand_ws.to_excel(pb_file, sheet_name="Random Duplicates")
        except ValueError:
            pass
        # Other Duplicates
        try:
            other_ws = pd.DataFrame.from_dict(self.duplicated_other, orient='index')
            other_ws.to_excel(pb_file, sheet_name="Other Duplicates")
        except ValueError:
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
            org_ms_count = pd.DataFrame.from_dict(org_gene_ms_count, orient='index')
            org_ms_count.to_excel(pb_file, sheet_name="Missing Genes Count")
            org_ms = pd.DataFrame.from_dict(org_gene_ms, orient='index')
            org_ms.to_excel(pb_file, sheet_name="Missing Genes by Org")
        except ValueError:
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
            gene_ms_count = pd.DataFrame.from_dict(gene_org_ms_count, orient='index')
            gene_ms_count.to_excel(pb_file, sheet_name="Missing Organisms Count")
            gene_ms = pd.DataFrame.from_dict(gene_org_ms, orient='index')
            gene_ms.to_excel(pb_file, sheet_name="Missing Organisms by Genes")
        except ValueError:
            pass
        pb_file.save()



