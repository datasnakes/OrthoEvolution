"""Class for managing, downloading and extracting features from genbank files."""

import os
import shutil
from pathlib import Path
from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Datasnakes.Tools import LogIt
from Datasnakes.Orthologs.utils import attribute_config
from Datasnakes.Orthologs.Blast.blastn_comparative_genetics import CompGenBLASTn
from Datasnakes.Tools.utils.other_utils import makedirectory
from Datasnakes.Orthologs.Blast.comparative_genetics_objects import CompGenObjects


class GenBank(object):
    def __init__(self, project, project_path=None, solo=False, multi=True, archive=False, min_fasta=True, blast=CompGenBLASTn, **kwargs):
        """
        This class will handle GenBank files in various ways.  It allows for refseq-release .gbff files to be downloaded
        from NCBI and uploaded to a BioSQL database (biopython).  Single .gbk files can be downloaded from the .gbff,
        and uploaded to a custom BopSQL database for faster acquisition of GenBank data.

        :param project:  The name of the project.
        :param project_path: The relative path to the project.
        :param solo:  A flag for adding single fasta files.
        :param multi:  A flag for adding multi-fasta files.
        :param archive: A flag for archiving current GenBank Data.  # TODO
        :param min_fasta: A flag for minimizing FASTA file headers.
        :param blast:  The blast parameter is used for composing various
                       Orthologs.Blast classes.  Can be a class, a dict,
                       or none.
        :returns:  .gbff files/databases, .gbk files/databases, and FASTA files.
        """

        # TODO-ROB: Change the way the file systems work.
        self.project = project
        self.solo = solo
        self.multi = multi
        self.min_fasta = min_fasta
        self.genbanklog = LogIt().default(logname="GenBank", logfile=None)

        # Configuration of class attributes
        add_self = attribute_config(self, composer=blast, checker=CompGenBLASTn, checker2=CompGenObjects,
                                    project=project, project_path=project_path)
        for var, attr in add_self.__dict__.items():
            setattr(self, var, attr)

        # Configuration
        self.target_gbk_db_path = self.user_db / Path(self.project)
        Path.mkdir(self.target_gbk_db_path, parents=True, exist_ok=True)

        # Make a list of BioSQL database(.db) files that contain GenBank info
        self.db_files_list = []
        for FILE in os.listdir(str(self.ncbi_db_repo)):
            if FILE.endswith('.db'):
                self.db_files_list.append(str(FILE))

    @staticmethod
    def name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode):
        """
        Provide a uniquely named FASTA file:
        * Coding sequence:
            * Single - "<path>/<gene>_<organism><feat_type_rank>.<extension>"
            * Multi  - "<path>/<gene><feat_type_rank>.<extension>"
        * Other:
            * Single - "<path>/<gene>_<organism>_<feat_type_rank>.<extension>"
            * Multi  - "<path>/<gene>_<feat_type_rank>.<extension>"

        :param path:  The path where the file will be made.
        :param gene:  The gene name.
        :param org:  The organism name.
        :param feat_type:  The type of feature from the GenBank record.  (CDS, UTR, misc_feature, variation, etc.)
        :param feat_type_rank:  The feature type  + the rank.  (There can be multiple misc_features and variations)
        :param extension:  The file extension.  (".ffn", ".faa", ".fna", ".fasta")
        :param mode:  The mode ("w" or "a") for writing the file.  Write to a solo-FASTA file.  Append a multi-FASTA
                      file.
        :return:  The uniquely named FASTA file.
        """

        # Create path variables.  (typically raw_data/<gene>/GENBANK
        feat_path = path
        # Create a format-able string for file names
        if feat_type_rank is "CDS":
            single = '%s_%s%s%s'
            multi = '%s%s%s'
        else:
            single = '%s_%s_%s%s'
            multi = '%s_%s%s'
        # Create different names based on fasta file type
        if mode == 'w':
            file_path = feat_path / Path(single % (gene, org, feat_type_rank, extension))
        elif mode == 'a':
            file_path = feat_path / Path(multi % (gene, feat_type_rank, extension))

        # Make the base directory and return an open file.
        makedirectory(feat_path)
        file = open(file_path, mode)
        return file

    @staticmethod
    def protein_gi_fetch(feature):
        """
        Retrieve the protein gi number.

        :param feature:  Search the protein feature for the GI number.
        :return:  The protein GI number as a string.
        """

        # Find the protein gi number under the features qualifiers.
        for x in feature.qualifiers:
            if 'GI' in x:
                head, sup, p_gi = x.partition(':')
                return p_gi

    def create_post_blast_gbk_records(self, org_list, gene_dict):
        """
        After a blast has completed and the accession numbers have been compiled into an accession file, this class
        searches a local NCBI refseq release database composed of GenBank records.  This method will create a single
        GenBank file (.gbk) for each ortholog with an accession number.  The create_post_blast_gbk_records is only callable if the
        the instance is composed by one of the Blast classes.  This method also requires an NCBI refseq release
        database to be set up with the proper GenBank Flat Files (.gbff) files.

        :param org_list:  List of organisms
        :param gene_dict:  A nested dictionary for accessing accession numbers.
                           (e.g. gene_dict[GENE][ORGANISM} yields an accession number)
        :return:  Does not return an object, but it does create all the proper genbank files.
        """
        # Parse the tier_frame_dict to get the tier
        for G_KEY, G_value in self.tier_frame_dict.items():
            tier = G_KEY
            # Parse the tier based transformed dataframe to get the gene
            for GENE in self.tier_frame_dict[tier].T:
                # Parse the organism list to get the desired accession number
                for ORGANISM in org_list:
                    accession = str(gene_dict[GENE][ORGANISM])
                    accession, sup, version = accession.partition('.')
                    # When parsing a GenBank database, the version needs to be removed.
                    accession = accession.upper()
                    server_flag = False
                    # Search the databases and create a GenBank file.
                    self.get_gbk_file(accession, GENE, ORGANISM, server_flag=server_flag)

    def get_gbk_file(self, accession, gene, organism, server_flag=None):
        """
        This method searches a GenBank database for a target accession number.  Generally, there are multiple GenBank
        database files for one NCBI dataset.

        :param accession: Accession number of interest without the version.
        :param gene: Target gene of the accession number parameter.
        :param organism: Target organism of the accession number parameter.
        :return: This function searches through the given NCBI databases (created by uploading NCBI refseq .gbff files
                 to a BioPython BioSQL database) and creates single GenBank files.  This function can be used after a
                 blast or on its own.  If used on it's own then the NCBI .db files must be manually moved to the proper
                 directories.
        """

        gene_path = self.raw_data / Path(gene) / Path('GENBANK')
        Path.mkdir(gene_path, parents=True, exist_ok=True)

        # Parse each database to find the proper GenBank record
        for FILE in self.db_files_list:
            db_file_path = self.ncbi_db_repo / Path(FILE)
            # Stop searching if the GenBank record has been created.
            if server_flag is True:
                break
            server = BioSeqDatabase.open_database(driver='sqlite3', db=str(db_file_path))
            # Parse the sub-databases
            for SUB_DB_NAME in server.keys():
                db = server[SUB_DB_NAME]
                try:
                    record = db.lookup(accession=accession)
                    gbk_file = '%s_%s.gbk' % (gene , organism)
                    gbk_file_path = gene_path / Path(gbk_file)
                    with open(gbk_file_path, 'w') as GB_file:
                        GB_file.write(record.format('genbank'))
                        self.genbanklog.info(GB_file.name, 'created')
                    # TODO-ROB:  Add quality control method here
                    # self.gbk_quality_control()
                    # Stop searching if the GenBank record has been created.
                    server_flag = True
                    break
                except IndexError:
                    self.genbanklog.critical('Index Error in %s.  Moving to the next database...' % SUB_DB_NAME)
                    continue

        # If the file has not been created after searching, then raise an error
        if server_flag is not True:
            self.genbanklog.critical("The GenBank file was not created for %s (%s, %s)." % (accession, gene, organism))
            raise FileNotFoundError

    def gbk_quality_control(self, gbk_file, gene, organism):
        print(self)
        print("Add this method to GenBank/utils.py")
        print("Check the GenBank file for proper data")

    def gbk_upload(self):
        """
        This method is only usable after creating GenBank records with this class.  It uploads a BioSQL databases with
        target GenBank data (.gbk files).  This creates a compact set of data for each project.

        :return:  Does not return an object, but creates a database for each gene-tier in the dataset.
        """

        t_count = 0
        # Parse the tier dictionary
        for TIER in self.tier_frame_dict.keys():
            db_name = str(TIER) + '.db'
            db_file_path = self.target_gbk_db_path / Path(db_name)
            # Create the db file if it exists
            if os.path.isfile(str(db_file_path)) is False:
                self.genbanklog.warn('Copying Template BioSQL Database...  This may take a few minutes...')
                shutil.copy2('Template_BioSQL_DB.db', str(db_file_path))

            # If it already exists then the database is bad, or needs to be update.  Delete it.
            else:
                # TODO-ROB: This part is broken until the template db creation and management is added
                os.remove(str(db_file_path))
                self.genbanklog.warn('Copying Template BioSQL Database...  This may take a few minutes...')
                shutil.copy2('Template_BioSQL_DB.db', str(db_file_path))

            server = BioSeqDatabase.open_database(driver='sqlite3', db=str(db_file_path))
            gene_path = self.raw_data
            # Parse the raw_data folder to get the name of each gene.
            for GENE in os.listdir(str(gene_path)):
                sub_db_name = GENE
                genbank_path = gene_path / Path(GENE) / Path('GENBANK')
                # Parse the GenBank file names for each gene in order to upload them to a custom BioSQL database
                for FILE in os.listdir(str(genbank_path)):
                    # Try to load the database.
                    try:
                        if sub_db_name not in server.keys():
                            server.new_database(sub_db_name)
                        db = server[sub_db_name]
                        count = db.load(SeqIO.parse(FILE, 'genbank'))
                        server.commit()
                        self.genbanklog.info('Server Commited %s' % sub_db_name)
                        self.genbanklog.info('%s database loaded with %s.' % (db.dbid, FILE))
                        self.genbanklog.info("That file contains %s genbank records." % str(count))
                        t_count = t_count + count
                        self.genbanklog.info('The total number of files loaded so far is %i.' % t_count)
                    # If the database cannot be loaded then rollback the server and raise an error.
                    except BaseException:
                        server.rollback()
                        # Try to delete the sub database and commit
                        try:
                            del server[sub_db_name]
                            server.commit()
                        # If it cannot be deleted then raise an error.
                        except BaseException:
                            raise
                        raise

    def get_fasta_files(self, acc_dict, db=True):
        """
        This method creates FASTA files for every GenBank record in the accession number dictionary.  It can search
        through a BioSQL database or it can crawl a directory for .gbk files.

        :param acc_dict:  An accession dictionary like the one created by CompGenObjects.
        :param db:  A flag that determines whether or not to use the custom BioSQL database or to use .gbk files.
        :return:  Does not return an object, but it will return a set of FASTA files for each GenBank record.
        """

        # Get FASTA files from the BioSQL GenBank databases.
        if db is True:
            # Parse the directory that contains the databases for the project of interest.
            for database in os.listdir(str(self.target_gbk_db_path)):
                server = BioSeqDatabase.open_database(driver="sqlite3", db=database)
                try:
                    for db_name in server.keys():
                        db = server[db_name]
                        # For each GenBank record in the database write a set of FASTA files.
                        for item in db.keys():
                            record = db.lookup(item)
                            self.write_fasta_files(record, acc_dict)
                            self.genbanklog.info("FASTA files for %s created from BioSQL database." % item)
                except:
                    raise()
        # Get FASTA files from the GenBank files.
        # TODO-ROB change this.  Broken by new directory structure
        # TODO-ROB directory looks like /raw_data/Gene_1/GENBANK/*.gbk
        elif db is False:
            # Parse the directory that contain the GenBank records for the project of interest.
            for root, dirs, gbk_files in os.walk(str(self.target_gbk_files_path)):
                # For each genbank record write a set of FASTA files.
                for gbk_file in gbk_files:
                    if Path(gbk_file).suffix is '.gbk':
                        record = SeqIO.read(gbk_file, 'genbank')
                        self.write_fasta_files(record, acc_dict)
                        self.genbanklog.info("FASTA files for %s created." % gbk_file)

    def write_fasta_files(self, record, acc_dict):
        """
        This method initializes a FASTA file by creating a dictionary for formatting the FASTA header and the following
        sequence.

        :param record:  A GenBank record created by BioPython.
        :param acc_dict:  The accession dictionary from the CompGenObjects class.
        :return:
        """

        feat_type_list = []
        for feature in record.features:
            # ############ Set up variables to use for dictionary values ############# #
            # Basic variables.
            accession = record.id
            gene = acc_dict[accession][0]
            organism = acc_dict[accession][1]
            # Variable for minimalistic FASTA files.
            genus, sep, species = organism.partition('_')
            min_org = str(''.join([genus[0], sep, species[0:28]]))

            # Keep a list of feature types to identify duplicates (for naming the FASTA files).
            # The first iteration of the feature type contains no number.
            # The following iterations are concatenated with a number.
            feat_type = str(feature.type)
            feat_type_list.append(feat_type)
            duplicate_num = feat_type_list.count(feat_type)
            if duplicate_num == 1:
                feat_type_rank = feat_type
            else:
                feat_type_rank = feat_type + str(duplicate_num)
            # ############ End ############# #

            # ######### Create a dictionary and format FASTA file entries. ######### #
            fmt = {
                'na_gi': str(record.annotations['gi']),
                'aa_gi': str(self.protein_gi_fetch(feature)),
                'na_acc_n': str(accession),
                'aa_acc_n': str(feature.qualifiers['protein_id'][0]),
                'na_description': str(record.description),
                'aa_description': str(feature.qualifiers['product'][0]),
                'na_seq': str(feature.extract(record.seq)),
                'aa_seq': str(feature.qualifiers['translation'][0]),
                'na_misc_feat': str(feature.qualifiers['note'][0]),
                'org': str(organism),
                'gene': str(gene),
                'min_org': str(min_org),
                'feat_type': str(feat_type),
                'feat_type_rank': str(feat_type_rank),
                'path': str(self.raw_data / Path(gene) / Path('GENBANK'))
            }
            # Set up minimalistic FASTA headers and sequence entries for Nucleic Acid and Amino Acid sequences.
            na_entry = ">{min_org}\n{na_seq}\n".format(**fmt)
            aa_entry = ">{min_org}\n{aa_seq}\n".format(**fmt)
            # For full FASTA headers/sequences set min_fasta to False
            if self.min_fasta is False:
                na_entry = ">gi|{na_gi}|ref|{na_acc_n}| {na_description}\n{na_seq}\n".format(**fmt)
                aa_entry = ">gi|{aa_gi}|reg|{aa_acc_n}| {aa_description} {org}\n{aa_seq}\n".format(**fmt)
            # ######### End ######### #

            # ############ Write desired FASTA files ############ #
            if self.solo is True:
                self.solo_fasta(na_entry, aa_entry, fmt)
            if self.multi is True:
                self.multi_fasta(na_entry, aa_entry, fmt)

    def solo_fasta(self, na_entry, aa_entry, fmt):
        """
        This method writes a sequence of a feature to a uniquely named file using a dictionary for formatting.

        :param na_entry:  A string representing the Nucleic Acid sequence data in FASTA format.
        :param aa_entry:  A string representing the Amino Acid sequence data in FASTA format.
        :param fmt:  A dictionary for formatting the FASTA entries and the file names.
        :return:  Does not return an object, but creates single entry FASTA files.
        """
        mode = 'w'

        # Create the desired variables from the formatter dictionary.
        feat_type = fmt['feat_type']
        feat_type_rank = fmt['feat_type_rank']
        path = fmt['path']
        gene = fmt['gene']
        org = fmt['org']

        if feat_type == "CDS":
            # Create a .ffn file (FASTA for Coding Nucleic Acids)
            extension = '.ffn'
            file = self.name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()
            # Create a .faa file (FASTA for Amino Acids)
            extension = '.faa'
            file = self.name_fasta_file(path, gene, org, 'Protein', feat_type_rank, extension, mode)
            file.write(aa_entry)
            file.close()

        elif feat_type == "misc_feature":
            # Create a custom entry for miscellaneous features.
            na_entry = ">gi|{na_gi}|ref|{na_acc_n}| {na_description} Feature: {na_misc_feat}\n{na_seq}\n".format(**fmt)
            # Creates .fna files (generic FASTA file for Nucleic Acids)
            extension = '.fna'
            file = self.name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()

        elif feat_type != "variation":
            # Creates .fasta files (generic FASTA file)
            extension = '.fasta'
            file = self.name_fasta_file(path, gene, org, 'Other', feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()

    def multi_fasta(self, na_entry, aa_entry, fmt):
        """
        This method appends an othologous sequence of a feature to a uniquely named file using a dictionary for
        formatting.

        :param na_entry:  A string representing the Nucleic Acid sequence data in FASTA format.
        :param aa_entry:  A string representing the Amino Acid sequence data in FASTA format.
        :param fmt:  A dictionary for formatting the FASTA entries and the file names.
        :return:  Does not return an object, but creates or appends to a multi entry FASTA file.
        """
        mode = 'a'

        # Create the desired variables from the formatter dictionary.
        feat_type = fmt['feat_type']
        feat_type_rank = fmt['feat_type_rank']
        path = fmt['path']
        gene = fmt['gene']
        org = fmt['org']

        if feat_type == "CDS":
            # Create a MASTER .ffn file (multi-FASTA file for Coding Nucleic Acids)
            extension = '.ffn'
            file = self.name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()
            # Create a MASTER .faa file (multi-FASTA file for Amino Acids)
            extension = '.faa'
            file = self.name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode)
            file.write(aa_entry)
            file.close()
        elif feat_type == "misc_feature":
            na_entry = ">gi|{na_gi}|ref|{na_acc_n}| {na_description} Feature: {na_misc_feat}\n{na_seq}\n".format(**fmt)
            # Creates .fna files (generic FASTA file for Nucleic Acids)
            extension = '.fna'
            file = self.name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode)
            file.write(na_entry)
            file.close()
