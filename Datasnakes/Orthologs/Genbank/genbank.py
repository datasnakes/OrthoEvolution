"""Class for managing, downloading and extracting features from genbank files."""

import os
import shutil
from pathlib import Path
from BioSQL import BioSeqDatabase
from Bio import SeqIO
from Datasnakes.Orthologs.Blast.blastn import BLASTn

# TODO-ROB:  REMOVED Tier Based Directory System.  Only add tier directories at the end of analysis in the users data folder

class GenBank(BLASTn):
    """Class for managing, downloading and extracting features from genbank files."""

    def __init__(self, solo=False, multi=True, archive=False, min_fasta=True, **kwargs):
        """Handle genbank files in various ways for the Orthologs Project.

        :param ncbi_db_repo: A path to the .db files of interest.  These
        .db files were created by downloading NCBI's refseq-GenBank
        flat files, and then uploading these files to a .db file
        that uses biopython's BioSQL database schema.
        :param gbk_path: The path used for the custom .db files and
        the single entry GenBank files.
        use:
            Path(gbk_path) / Path(target_gbk_files_path);
            Path(gbk_path) / Path(db_files)

        :param solo: Boolean for adding single fasta files.
        :param multi:  Boolean for adding multi-fasta files.
        """

        # TODO-ROB: Change the way the file systems work.
        super().__init__(**kwargs)
        self.gbk_path = self.raw_data / Path('GENBANK')
        self.target_gbk_files_path = self.gbk_path / Path('gbk')  # Deprecate this
        self.target_gbk_files_path.mkdir(parents=True, exist_ok=True)

        self.target_fasta_files = self.gbk_path / Path('fasta')  # Deprecated this
        self.target_fasta_files.mkdir(parents=True, exist_ok=True)

        self.target_gbk_db_path = self.user_db / Path(self.project)  # Configuration.  Create configuration script.
        # TODO-ROB: Configure GenBank function
        self.target_gbk_db_path.mkdir(parents=True, exist_ok=True)
        self.solo = solo
        self.multi = multi
        self.min_fasta = min_fasta

    @staticmethod
    def name_fasta_file(path, gene, org, feat_type, feat_type_rank, extension, mode):
        """Provide a unique name for the fasta file."""

        # Create path variables.
        feat_path = path / Path(feat_type) / Path(gene)
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
        feat_path.mkdir(parents=True, exist_ok=True)
        file = open(file_path, mode)
        return file

    @staticmethod
    def protein_gi_fetch(feature):
        """Retrieve the protein gi number."""

        # Find the protein gi number under the features qualifiers.
        for x in feature.qualifiers:
            if 'GI' in x:
                head, sup, p_gi = x.partition(':')
                return p_gi

    def get_gbk_files(self, org_list, gene_dict):
        """
        :param org_list:  List of organisms
        :param gene_dict:  A nested dictionary for accessing accession numbers.
        (e.g. gene_dict[GENE][ORGANISM} yields an accession number)
        :return:  Does not return an object, but it does create all the proper
        genbank files.
        """
        # Make a list of BioSQL database(.db) files that contain GenBank info
        db_files_list = []
        for FILE in os.listdir(str(self.ncbi_db_repo)):
            if FILE.endswith('.db'):
                db_files_list.append(str(FILE))

        # Parse the tier_frame_dict to get the tier
        for G_KEY, G_value in self.tier_frame_dict.items():
            tier = G_KEY
            # Parse the tier based transformed dataframe to get the gene
            for GENE in self.tier_frame_dict[tier].T:
                gene_path = self.target_gbk_files_path / Path(GENE)
                Path.mkdir(gene_path)
                # Parse the organism list to get the desired accession number
                for ORGANISM in org_list:
                    accession = str(gene_dict[GENE][ORGANISM])
                    accession, sup, version = accession.partition('.')
                    accession = accession.upper()
                    server_flag = False
                    # Parse each database to find the proper GenBank record
                    for FILE in db_files_list:
                        db_file_path = self.ncbi_db_repo / Path(FILE)
                        if server_flag is True:
                            break
                        server = BioSeqDatabase.open_database(driver='sqlite3', db=str(db_file_path))
                        # Parse the sub-databases
                        for SUB_DB_NAME in server.keys():
                            db = server[SUB_DB_NAME]
                            try:
                                record = db.lookup(accession=accession)
                                gbk_file = '%s_%s.gbk' % (GENE, ORGANISM)
                                gbk_file_path = gene_path / Path(gbk_file)
                                with open(gbk_file_path, 'w') as GB_file:
                                    GB_file.write(record.format('genbank'))
                                    print(GB_file.name, 'created')
                                server_flag = True
                                break
                            except IndexError:
                                print('Index Error in %s.  Moving to the next database...' % SUB_DB_NAME)
                                continue

    def get_fasta_files(self, acc_dict, db=True):
        """Create FASTA files for every GenBank record."""

        # Get FASTA files from the BioSQL GenBank databases.
        if db is True:
            for database in os.listdir(str(self.target_gbk_db_path)):
                server = BioSeqDatabase.open_database(driver="sqlite3", db=database)
                try:
                    for db_name in server.keys():
                        db = server[db_name]
                        for item in db.keys():
                            record = db.lookup(item)
                            self.write_fasta_files(record, acc_dict)
                except:
                    raise()
        # Get FASTA files from the GenBank files.
        elif db is False:
            for root, dirs, gbk_files in os.walk(str(self.target_gbk_files_path)):
                for gbk_file in gbk_files:
                    if Path(gbk_file).suffix is '.gbk':
                        record = SeqIO.read(gbk_file, 'genbank')
                        self.write_fasta_files(record, acc_dict)

    def gbk_upload(self):
        """Upload the BioSQL database with genbank data."""
        t_count = 0
        Path.mkdir(self.target_gbk_db_path)
        for TIER in self.tier_frame_dict.keys():
            db_name = str(TIER) + '.db'
            db_file_path = self.target_gbk_db_path / Path(db_name)
            if os.path.isfile(str(db_file_path)) is False:
                print('Copying Template BioSQL Database...  This may take a few minutes...')
                # TODO-ROB:  Create a utility function for creating BioSQL databases
                shutil.copy2('Template_BioSQL_DB.db', str(db_file_path))
            else:
                # This part is broken until the template db creation and management is added
                os.remove(str(db_file_path))
                print('Copying Template BioSQL Database...  This may take a few minutes...')
                shutil.copy2('Template_BioSQL_DB.db', str(db_file_path))

            server = BioSeqDatabase.open_database(driver='sqlite3', db=str(db_file_path))
            gene_path = self.target_gbk_files_path
            for GENE in os.listdir(str(gene_path)):
                sub_db_name = GENE
                gene_path = gene_path / Path(GENE)
                for FILE in os.listdir(str(gene_path)):
                    try:
                        if sub_db_name not in server.keys():
                            server.new_database(sub_db_name)
                        db = server[sub_db_name]
                        count = db.load(SeqIO.parse(FILE, 'genbank'))
                        server.commit()
                        print('Server Commited %s' % sub_db_name)
                        print('%s database loaded with %s.' % (db.dbid, FILE))
                        print("That file contains %s genbank records." % str(count))
                        t_count = t_count + count
                        print('The total number of files loaded so far is %i.' % t_count)
                    except BaseException:
                        server.rollback()
                        try:
                            del server[sub_db_name]
                            server.commit()
                        except BaseException:
                            raise
                        raise

    def write_fasta_files(self, record, acc_dict):
        """Initialize writing a fasta sequence or feature to a file."""
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
                'path': self.target_fasta_files
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
        """Write a feature to a file."""
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
        """Write multiple fasta files."""
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
