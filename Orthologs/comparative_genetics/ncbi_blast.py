import time
import os
from pathlib import Path
import pandas as pd
from Orthologs.comparative_genetics.comp_gen import CompGenAnalysis as CGA
# from Orthologs.manager.logit.logit import LogIt


class BLASTAnalysis(CGA):
    def __init__(self, repo, user, project, research, research_type,
                 template=None, taxon_file=None, post_blast=False):
        """Inherit from the CompGenAnalysis class.  If the BLAST was cut short,
        then a build_file is to be used."""
        super().__init__(repo=repo, user=user, project=project, research=research, research_type=research_type,
                         acc_file=template, taxon_file=taxon_file, post_blast=post_blast, save_data=True, hgnc=False)
        # TODO-ROB: Inherit or add variable for logger class
        # TODO-ROB Add Mana directories
        # Private variables
        self.__home = os.getcwd()
        # name, ext = os.path.splitext(template)
        # new_name = "%s_template%s" % (self.project, ext)
        # os.rename(template, new_name)
        self.template_filename = self.acc_filename
        self.template_path = self.data / Path(self.template_filename)
        self.__taxon_filename = taxon_file
        self.__post_blast = post_blast
        self.__save_data = True

        # Initialize Logging
        df = LogIt()
        self.blastn_log = df.blastn()
        self.postblast_log = df.post_blast()
        self.config_log = df.config(self.user_log / Path('BLAST.log'))
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

        # Initialize a data frame to add accession files to
        self.building = self.raw_acc_data
        del self.building['Tier']
        del self.building['Homo_sapiens']
        self.building = self.building.set_index('Gene')
        build_file = str(template[:-4] + 'building.csv')
        self.__building_filename = build_file
        self.__building_file_path = self.data / Path(self.__building_filename)

        # Initialize some timing files to add times to
        self.building_time = self.building
        time_file = build_file.replace('building.csv', 'building_time.csv')
        self.__building_time_filename = time_file
        self.__building_time_file_path = self.data / Path(self.__building_time_filename)

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
                self.blastn_log.critical("Existing Accession Organism: %s" % self.acc_dict[existing][0])
                self.blastn_log.critical("Existing Accession Organism: %s" % self.acc_dict[existing][1])
                self.blastn_log.critical("***********************************************************************")
                raise ValueError("The queried ACCESSION (%s) does not match the existing ACCESSION (%s).  Please see"
                                 "the log file." % (accession, existing))

        self.building.set_value(gene, organism, accession)
        temp = self.building.reset_index()
        temp.insert(0, 'Tier', self.df['Tier'])
        if self.__save_data is True:
            temp.to_csv(self.__building_file_path)

    def add_blast_time(self, gene, organism, start, end):
        """Takes the start/end times and adds in to the building_time 
        dataframe, and also adds to the csv file."""
        elapsed_time = end - start
        # Edit the data frame
        self.building_time.set_value(gene, organism, elapsed_time)
        temp = self.building_time.reset_index()
        temp.insert(0, 'Tier', self.df['Tier'])
        if self.__save_data is True:
            temp.to_csv(self.__building_time_file_path)

    def post_blast_analysis(self, acc_file):
        # TODO-ROB Make this into the CGA.make_excel_files() method
        accession_data = CGA(acc_file=acc_file, post_blast=True)
        self.postblast_log.info('*************************POST BLAST ANALYSIS START*************************\n\n\n')

        missing_gene = accession_data.missing_dict['genes']
        missing_orgs = accession_data.missing_dict['organisms']
        orgs = accession_data.org_list

        # Create and write dictionaries & dataframes to excel file
        if missing_gene['count'] <= 0 and missing_orgs['count'] <= 0:
            # Log that the blast had full coverage
            self.postblast_log.info('There are no missing accession numbers for any gene or organism.')
            self.postblast_log.info('Post blastn analysis is complete.')
        else:
            self.post_blast_log.info('There are missing accessions. This data will be written an excel file.')

            # Set up the excel file
            excel_file = pd.ExcelWriter(processed + 'karg_missing_genes_data.xlsx')

            # This is the data frame for the names of missing genes by organism
            frame = pd.DataFrame.from_dict(no_acc, orient='index', dtype=str)
            frame = frame.transpose()
            frame.to_excel(excel_file, sheet_name="Missing Genes by Org", index=False)

            # This is the data frame for the number of missing genes by organism
            frame2 = pd.DataFrame.from_dict(num_missing, orient='index')
            frame2.to_excel(excel_file, sheet_name="# of Genes missing by Org", header=False)

            # This is the data frame for the number of missing ORGANISMS by gene
            frame3 = pd.DataFrame(miss_per_gene)
            combine = [original.Tier, frame3]  # Combine the tier column and data column
            frame3 = pd.concat(combine, axis=1) # Concatenate the dataframes
            frame3.to_excel(excel_file, sheet_name="# of Orgs missing by Gene", header=False)

            # Save the excel file
            excel_file.save()
            post_blast_log.info('Your file, karg_missing_genes_data.xlsx, has been created and saved.')

            post_blast_log.info('Post blastn analysis is complete.')


