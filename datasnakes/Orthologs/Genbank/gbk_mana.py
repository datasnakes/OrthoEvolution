"""
This class is designed to aid in extracting genbank files from a prefilled
database.
"""
import os
import shutil
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from datasnakes.Manager.utils.mana import UserMana as UM
from datasnakes.Orthologs.CompGenetics import CompGenAnalysis as CGA

# TODO-ROB Add a progress bar to the pipeline
# TODO-ROB make code more versatile for multiple projects or even single
# queries


class Database2Genbank(CGA, UM):
    """
    Extract target GenBank files from the database files that were created using
    ftp2db.
    """

    def __init__(self, repo, user, project, m_file, new_db=True):
        # if new db is False then update the genbank files
        self.GenBank_update = new_db
        # TODO-ROB Update the CGA.  I have yet to write functions for saving
        # the data post blast
        CGA.__init__(repo=repo, user=user, project=project, acc_file=m_file, save_data=False)
        UM.__init__(repo=repo, user=user, porject=project, database=['genbank'], new_db=new_db)

    def gbk_gather(self):
        """Extract/download the genbank files from the database.
        """
        os.chdir(self.databases)
        db_name = []

        for file in os.listdir(os.getcwd()):
            if file.endswith('.db'):
                db_name.append(str(file))

        for G_key, G_value in self.what.tier_frame_dict.items():
            Tier = G_key
            os.chdir(self.path)
            os.mkdir(Tier)
            os.chdir(Tier)
            Tier_path = os.getcwd()
            for Gene in self.what.tier_frame_dict[Tier].T:
                os.chdir(Tier_path)
                os.mkdir(Gene)
                os.chdir(Gene)
                for Organism in self.what.org_list:
                    Accession = str(self.what.gene_dict[Gene][Organism])
                    Accession, Sup, Version = Accession.partition('.')
                    Accession = Accession.upper()
                    server_flag = False
                    for name in db_name:
                        if server_flag is True:
                            break
                        name = str(name)
                        server = BioSeqDatabase.open_database(
                            driver='sqlite3', db=where.VERT_MAM + ('/Databases/%s' % name))
                        for sub_db_name in server.keys():
                            db = server[sub_db_name]

                            try:
                                record = db.lookup(accession=Accession)
                                with open('%s_%s.gbk' % (Gene, Organism), 'w') as GB_file:
                                    GB_file.write(record.format('genbank'))
                                    print(GB_file.name, 'created')
                                server_flag = True
                                break
                            except IndexError:
                                print('Index Error')
                                continue

    def gbk_upload(self):
        """
        Upload the BioSQL database with genbank data.
        """
        t_count = 0
        os.chdir(self.path)
        print(os.getcwd())
        if os.path.isdir(self.path + '/Databases') is False:
            os.mkdir('Databases')
        for tier in os.listdir(os.getcwd()):
            if tier == 'Databases':
                continue
            db_name = str(tier) + '.db'
            if os.path.isfile(self.path + '/Databases/' + db_name) is False:
                print('Copying Template BioSQL Database...  '
                      'This may take a few minutes...')
                shutil.copy2(where.Templates + '/Template_BioSQL_DB.db',
                             self.path + '/Databases/%s' % db_name)
            else:
                os.remove(self.path + '/Databases/' + db_name)
                print('Copying Template BioSQL Database...  '
                      'This may take a few minutes...')
                shutil.copy2(where.Templates + '/Template_BioSQL_DB.db',
                             self.path + '/Databases/%s' % db_name)

            server = BioSeqDatabase.open_database(
                driver='sqlite3', db=(
                    self.path + '/Databases/' + db_name))
            os.chdir(tier)
            for gene in os.listdir(os.getcwd()):
                os.chdir(gene)
                sub_db_name = gene
                for file in os.listdir(os.getcwd()):
                    try:
                        if sub_db_name not in server.keys():
                            server.new_database(sub_db_name)
                        db = server[sub_db_name]
                        count = db.load(SeqIO.parse(file, 'genbank'))
                        server.commit()
                        print('Server Commited %s' % sub_db_name)
                        print('%s database loaded with %s.' % (db.dbid, file))
                        print(
                            "That file contains %s genbank records." %
                            str(count))
                        t_count = t_count + count
                        print(
                            'The total number of files loaded so far is %i.' %
                            t_count)
                    except BaseException:
                        server.rollback()
                        try:
                            del server[sub_db_name]
                            server.commit()
                        except BaseException:
                            raise
                        raise
                os.chdir('..')
            os.chdir('..')

    def db2gbk_mkdir(self, path, p_list, update):
        """
        This simple function aids in the handling of directory creation.
        """
        if update is True:
            path = where.dir_archive(path, p_list)
        else:
            path = where.dir_make(path, p_list)
        return path
