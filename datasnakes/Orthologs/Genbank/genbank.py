import os
import shutil
from pathlib import Path
from BioSQL import BioSeqDatabase
from Bio import SeqIO


class GenBank(object):

    def __init__(self, ncbi_db_path, gbk_path, solo=False, multi=True, archive=False, min_fasta=True):
        """
        Base GenBank class that handles GenBank files in various ways
        for the Orthologs Project.
        :param ncbi_db_path: A path to the .db files of interest.  These
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
        self.ncbi_db_path = Path(ncbi_db_path)
        self.gbk_path = Path(gbk_path) / Path('genbank')
        self.target_gbk_files_path = self.gbk_path / Path('gbk')
        self.target_gbk_db_path = self.gbk_path / Path('db')
        self.target_fasta_files = self.gbk_path / Path('solo')
        self.solo = solo
        self.multi = multi
        self.min_fasta = min_fasta

    def get_gbk_files(self, tier_frame_dict, org_list, gene_dict):
        """Extract/download the genbank files from the database.
        """
        # Make a list of BioSQL database(.db) files that contain GenBank info
        db_files_list = []
        for FILE in os.listdir(str(self.ncbi_db_path)):
            if FILE.endswith('.db'):
                db_files_list.append(str(FILE))

        for G_KEY, G_value in tier_frame_dict.items():
            tier = G_KEY
            tier_path = self.target_gbk_files_path / Path(tier)
            Path.mkdir(tier_path, parents=True, exist_ok=True)
            for GENE in tier_frame_dict[tier].T:
                gene_path = tier_path / Path(GENE)
                Path.mkdir(gene_path)
                for ORGANISM in org_list:
                    accession = str(gene_dict[GENE][ORGANISM])
                    accession, sup, version = accession.partition('.')
                    accession = accession.upper()
                    server_flag = False
                    for FILE in db_files_list:
                        db_file_path = self.ncbi_db_path / Path(FILE)
                        if server_flag is True:
                            break
                        server = BioSeqDatabase.open_database(driver='sqlite3', db=str(db_file_path))
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

    def gbk_upload(self):
        """
        Upload the BioSQL database with genbank data.
        """
        t_count = 0
        Path.mkdir(self.target_gbk_db_path)
        for TIER in os.listdir(str(self.target_gbk_files_path)):
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
            tier_path = self.target_gbk_files_path / Path(TIER)
            for GENE in os.listdir(str(tier_path)):
                sub_db_name = GENE
                gene_path = tier_path / Path(GENE)
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

    def get_fasta_files(self, acc_dict):
        """
        Create FASTA files for every GenBank record in the databases.
        """
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

    def write_fasta_files(self, record, acc_dict):
        min = self.min_fasta
        if self.solo is True:
            self.sol_fasta(record, acc_dict, min)
        if self.multi is True:
            self.multi_fasta(record, acc_dict, min)

        feat_type_list = []
        for feature in record.features:
            # TODO-ROB logic for self.multi
            m = 'multi'
            accession = record.id
            try:
                organism = acc_dict[accession][1]
            except KeyError:
                accession = str(record.id).lower()
                organism = acc_dict[accession][1]
            gi = record.annotations['gi']

            feat_type = str(feature.type)
            feat_type_list.append(feat_type)

            # Add this in the name_fasta_file definition
            # If there are duplicates of feature.type
            # then add a number to the end of the name
            x = feat_type_list.count(feat_type)  # Number of duplicates
            feat_type_rank = feat_type + str(x)  # Adding the number to the name of the feature.type
            if x == 1:
                feat_type_rank = feat_type
            ###########################

    def solo_fasta(self, record, acc_dict, min):
        print('')

        feat_type_list = []
        for feature in record.features:
            # TODO-ROB logic for self.multi
            m = 'multi'
            accession = record.id
            try:
                organism = acc_dict[accession][1]
            except KeyError:
                accession = str(record.id).lower()
                organism = acc_dict[accession][1]
            gi = record.annotations['gi']

            feat_type = str(feature.type)
            feat_type_list.append(feat_type)

            # Add this in the name_fasta_file definition
            # If there are duplicates of feature.type
            # then add a number to the end of the name
            x = feat_type_list.count(feat_type)  # Number of duplicates
            feat_type_rank = feat_type + str(x)  # Adding the number to the name of the feature.type
            if x == 1:
                feat_type_rank = feat_type
            if feat_type == "CDS":

                # Create a .ffn file (FASTA for Coding Nucleic Acids)
                self.file = self.name_fasta_file(accession, feat_type, feat_type_rank, extension='.ffn')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)))
                self.file.close()

                # Create a .faa file (FASTA for Amino Acids)
                self.file = self.name_fasta_file(accession, 'Protein', feat_type_rank, extension='.faa')
                gi_p = self.protein_gi_fetch(feature)
                self.file.write(">gi|" + str(gi_p) + "|ref|"
                                + str(feature.qualifiers['protein_id'][0]) + '| '
                                + str(feature.qualifiers['product'][0]) + ' ' + str(organism) + '\n'
                                + str(feature.qualifiers['translation'][0]))
                self.file.close()

            elif feat_type == "misc_feature":

                # Creates .fna files (generic FASTA file for Nucleic Acids)
                self.file = self.name_fasta_file(accession, feat_type, feat_type_rank, extension='.fna')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description
                                + ' Feature: ' + str(feature.qualifiers['note'][0]) + '\n'
                                + str(feature.extract(record.seq)))
                self.file.close()

            elif feat_type != "variation":
                print(feat_type)
                # Creates .fasta files (generic FASTA file)
                self.file = self.name_fasta_file(accession, 'Other', feat_type_rank, extension='.fasta')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

                self.file = self.name_fasta_file(accession, 'Other', feat_type_rank, extension='.fasta')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

    def multi_fasta(self, record, acc_dict, min):
        print('')

        feat_type_list = []
        for feature in record.features:
            # TODO-ROB logic for self.multi
            m = 'multi'
            accession = record.id
            try:
                organism = acc_dict[accession][1]
            except KeyError:
                accession = str(record.id).lower()
                organism = acc_dict[accession][1]
            gi = record.annotations['gi']

            feat_type = str(feature.type)
            feat_type_list.append(feat_type)

            # Add this in the name_fasta_file definition
            # If there are duplicates of feature.type
            # then add a number to the end of the name
            x = feat_type_list.count(feat_type)  # Number of duplicates
            feat_type_rank = feat_type + str(x)  # Adding the number to the name of the feature.type
            if x == 1:
                feat_type_rank = feat_type
            if feat_type == "CDS":

                # Create a MASTER .ffn file (multi-FASTA file for Coding Nucleic Acids)
                self.file = self.name_fasta_file(accession, m, feat_type_rank, extension='.ffn')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

                # Create a MASTER .faa file (multi-FASTA file for Amino Acids)
                self.file = self.name_fasta_file(accession, m, feat_type_rank, extension='.faa')
                gi_p = self.protein_gi_fetch(feature)
                self.file.write(">gi|" + str(gi_p) + "|ref|"
                                + str(feature.qualifiers['protein_id'][0]) + '| '
                                + str(feature.qualifiers['product'][0]) + ' ' + organism + '\n'
                                + str(feature.qualifiers['translation'][0]) + '\n')
                self.file.close()

            elif feat_type == "misc_feature":

                # Creates .fna files (generic FASTA file for Nucleic Acids)
                self.file = self.name_fasta_file(accession, m, feat_type_rank, extension='.fna')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description
                                + ' Feature: ' + str(feature.qualifiers['note'][0]) + '\n'
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

            elif feat_type != "variation":
                # TODO-ROB: fix for multifasta
                print(feat_type)
                # Creates .fasta files (generic FASTA file)
                self.file = self.name_fasta_file(accession, 'Other', feat_type_rank, extension='.fasta')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

                self.file = self.name_fasta_file(accession, 'Other', feat_type_rank, extension='.fasta')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

    def name_fasta_file(self, accession, feat_type, feat_type_rank, extension, acc_dict):
        gene = acc_dict[accession][0]
        organism = acc_dict[accession][1]
        # Create Handles for directories and file paths
        gene_feat_type_path = self.target_fasta_files / Path(feat_type) / Path(gene)
        file_path = gene_feat_type_path / Path('%s_%s_%s%s' % (gene, organism, feat_type_rank, extension))
        multi_file_path = gene_feat_type_path / Path('%s_%s%s' % (gene, feat_type_rank, extension))
        # Make the base directory
        gene_feat_type_path.mkdir(parents=True, exist_ok=True)

        if feat_type == 'MASTER':
            mode = 'a'
            self.file = open(multi_file_path, mode)
        else:
            mode = 'w'
            self.file = open(file_path, mode)
        return self.file

    def protein_gi_fetch(self, feature):
        for x in feature.qualifiers:
            if 'GI' in x:
                head, sup, gi_p = x.partition(':')
                return gi_p