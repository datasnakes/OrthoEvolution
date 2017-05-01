from manager.lister import *
from manager.logit.logit import LogIt
import time
import logging as log
from _datetime import datetime as d


class BLASTingTemplate(Lister):
    # TODO-ROB:  Look into subclasses and inheritance.  Which of the parameter below is necessary???
    def __init__(self, template=None, build_file=None, taxon_file=None, post_blast=False):
        """Inherit from the Lister class.  If the BLAST was cut short,
        then a build_file is to be used."""
        super().__init__(acc_file=template, taxon_file=taxon_file,
                         post_blast=post_blast, save_data=True, hgnc=False)
        # TODO-ROB: Inherit or add variable for logger class

        # Private variables
        self.__home = os.getcwd()
        self.__template_filename = template
        self.__acc_filename = build_file
        self.__taxon_filename = taxon_file
        self.__post_blast = post_blast
        self.__save_data = True

        # Initialize Logging
        df = LogIt()
        self.blastn_log = df.blastn()
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
        self.building = self.raw_data
        del self.building['Tier']
        del self.building['Homo_sapiens']
        self.building = self.building.set_index('Gene')
        self.building_time = self.building
        if build_file is None:
            build_file = str(template[:-4] + 'building.csv')
        time_file = build_file.replace('build.csv', 'building_time.csv') # TODO-ROB:  ADD time to the file name
        self.__building_filename = build_file
        self.__building_time_filename = time_file
        self.__building_file_path = Path(home) / Path('index') / Path(self.__building_filename)
        self.__building_time_file_path = Path(home) / Path('index') / Path(self.__building_time_filename)

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





