import os
from pathlib import Path
from pathlib import PurePath
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from dir_mana import dir_mana
from lister import Lister

home = os.getcwd()
project = "GPCR-Orthologs-Project"
user = "rgilmore"
where = dir_mana(home, project)
# Use lister() class here so that we can easily access our Master RNA Accession File
what = Lister('MAFV3.1.csv')  # Always make sure this file name is correct

class Gbk2Fasta(object):
    """
    This class accepts .gbk files or .db files (loaded with target .gbk files
    with BioPython schema).  The output is a group of FASTA files. The
    different types of regions that are specified by the GenBank files get
    their own directories.
    """

    def __init__(self, path=Path(os.getcwd()), db_flag=False, gbk_flag=False):
        # Where are the files?
        path = Path('C:\\Users\\rgilmore\\PycharmProjects\\GPCR-Orthologs-Project\\CODE\\1_Databases\\Vallender_Data')
        self.path = path

        # Parse .gbk files, .db files, or both
        os.chdir(self.path)
        for f_db in os.listdir():
            if '.db' in PurePath(f_db).suffix:
                self.db(f_db)
            elif '.gbk' in PurePath(f_db).suffix:
                self.gbk(f_db)

    def db(self, database):
        """
        Create FASTA files for every GenBank record in the database.
        """
        server = BioSeqDatabase.open_database(driver="sqlite3", db=database)
        try:
            for db_name in server.keys():
                db = server[db_name]
                for item in db.keys():
                    record = db.lookup(item)

                    self.write_fasta_file(record)
        except:
            raise()

    def gbk(self, file):
        """
        Write a group of fasta files from a group of .gbk files in the directory.
        The specified self.path must be an empty directory with .gbk files
        """
        record = SeqIO.read(file, 'genbank')
        self.write_fasta_file(record)

    def write_fasta_file(self, record):

        feature_list = []

        for feature in record.features:
            m = 'MASTER'
            accession = record.id
            try:
                organism = what.acc_dict[accession][1]
            except KeyError:
                accession = str(record.id).lower()
                organism = what.acc_dict[accession][1]
            gi = record.annotations['gi']

            feat = str(feature.type)
            feature_list.append(feature.type)

            # Add this in the name_fasta_file definition
            # If there are duplicates of feature.type
            # then add a number to the end of the name
            if feature.type in feature_list:
                x = feature_list.count(feature.type)  # Number of duplicates
                feat += str(x)  # Adding the number to the name of the feature.type

            if feature.type == "CDS":

                # Create a .ffn file (FASTA for Coding Nucleic Acids)
                file = self.name_fasta_file(accession, feature.type, feat, extension='.ffn')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)))
                self.file.close()

                # Create a .faa file (FASTA for Amino Acids)
                self.file = self.name_fasta_file(accession, 'Protein', feat, extension='.faa')
                gi_p = self.protein_gi_fetch(feature)
                self.file.write(">gi|" + str(gi_p) + "|ref|"
                                + str(feature.qualifiers['protein_id'][0]) + '| '
                                + str(feature.qualifiers['product'][0]) + ' ' + str(organism) + '\n'
                                + str(feature.qualifiers['translation'][0]))
                self.file.close()

                # Create a MASTER .ffn file (multi-FASTA file for Coding Nucleic Acids)
                self.file = self.name_fasta_file(accession, m, feat, extension='.ffn')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

                # Create a MASTER .faa file (multi-FASTA file for Amino Acids)
                self.file = self.name_fasta_file(accession, m, feat, extension='.faa')
                gi_p = self.protein_gi_fetch(feature)
                self.file.write(">gi|" + str(gi_p) + "|ref|"
                                + str(feature.qualifiers['protein_id'][0]) + '| '
                                + str(feature.qualifiers['product'][0]) + ' ' + organism + '\n'
                                + str(feature.qualifiers['translation'][0]) + '\n')
                self.file.close()

            elif feature.type == "misc_feature":

                # Creates .fna files (generic FASTA file for Nucleic Acids)
                self.file = self.name_fasta_file(accession, feature.type, feat, extension='.fna')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description
                                + ' Feature: ' + str(feature.qualifiers['note'][0]) + '\n'
                                + str(feature.extract(record.seq)))
                self.file.close()

                self.file = self.name_fasta_file(accession, m, feat, extension='.fna')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description
                                + ' Feature: ' + str(feature.qualifiers['note'][0]) + '\n'
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

            elif feature.type != "variation":
                print(feature.type)
                # Creates .fasta files (generic FASTA file)
                self.file = self.name_fasta_file(accession, 'Other', feat, extension='.fasta')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

                self.file = self.name_fasta_file(accession, 'Other', feat, extension='.fasta')
                self.file.write(">gi|" + gi + "|ref|" + accession + '| ' + record.description + "\n"
                                + str(feature.extract(record.seq)) + "\n")
                self.file.close()

    def name_fasta_file(self, accession, f_type, feature_no, extension):
        gene = what.acc_dict[accession][0]
        organism = what.acc_dict[accession][1]
        if os.path.isdir(self.path / Path(f_type)) is False:
            os.mkdir(self.path / Path(f_type))
        if os.path.isdir(self.path / Path(f_type) / Path(gene)) is False:
            os.mkdir(self.path / Path(f_type) / Path(gene))
        if f_type == 'MASTER':
            mode = 'a'
            self.file = open(self.path / Path(f_type) / Path(gene) / Path('%s_%s%s' % (gene, feature_no, extension)), mode)
        else:
            mode = 'w'
            self.file = open(self.path / Path(f_type) / Path(gene) / Path('%s_%s_%s%s' % (gene, organism, feature_no, extension)), mode)
        return self.file

    def protein_gi_fetch(self, feature):
        for x in feature.qualifiers:
            if 'GI' in x:
                head, sup, gi_p = x.partition(':')
                return gi_p