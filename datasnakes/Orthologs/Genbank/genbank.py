import os
import shutil
from pathlib import Path
from BioSQL import BioSeqDatabase
from Bio import SeqIO

class GenBank(object):

    def __init__(self, ncbi_db_path, gbk_path, fasta=False, multi=True, archive=False):
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
            
        :param fasta: Boolean for adding single fasta files.
        :param multi:  Boolean for adding multi-fasta files.
        """
        self.ncbi_db_path = Path(ncbi_db_path)
        self.gbk_path = Path(gbk_path) / Path('genbank')
        self.target_gbk_files_path = self.gbk_path / Path('gbk')
        self.target_gbk_db_path = self.gbk_path / Path('db')
        self.fasta = fasta
        self.multi = multi

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
