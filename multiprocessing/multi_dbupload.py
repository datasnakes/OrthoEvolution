import json
import os
import platform
import time

from Bio import SeqIO
from BioSQL import BioSeqDatabase
from mpi4py import MPI

from dir_mana import dir_mana
from lister import Lister

# Get child process information
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
machine = platform.node()

# Initialize project based functions and variables
home = os.getcwd()
project = "GPCR-Orthologs-Project"
user = "rgilmore"
where = dir_mana(home, project)
what = Lister('MAFV3.1.csv')

# Initialize the files and variables used to upload the proper database files
log_file = ''
# THe download list contains 8 lists with file names  e.g. [[list1], [list2]....[list8]]
temp_file = 'temp_file.txt'   # A JSON file (keys: 'logfile'; 'downloadlist', 'db_name', 'key)
# (values: 'file'; list of lists;, 'refseq-release-vertebrate-mammalian', 'vertebrate-mammalian)
temp_file_flag = True
process_number = rank
home = home + '/CODE/1_Databases/NCBI_Data/refseq/release/multiprocessing'

# Try to open the temporary file.  The same file is trying to be accessed by 8 different jobs
# TODO-ROB Make 8 different files with 0-7 labeled
while temp_file_flag is True:
    try:
        with open(home + '/' + temp_file, mode='r', encoding='utf-8') as temp_file:
            temp_var = json.load(temp_file)
            temp_file_flag = False
            print('JSON file loaded from process %s' % process_number)
    except:
        temp_file_flag = True
        print('Trying to open %s again from process %s' % ((home+'/'+temp_file), process_number))
        time.sleep(2)

# Initialize the dictionary for the download list.  'temp_var' already has the JSON data with mentioned key/value pairs
temp_var['db_name'] = str(temp_var['db_name'] + str(process_number))  # Each process get's it's own unique database
temp_var['small_list'] = temp_var['download_list'][rank-1]   # Each process gets it's own unique list of files
temp_var['log_file_rank'] = str(temp_var['log_file']) + str(rank)  # Each process gets it's own unique log file
ser_loc = where.DB
loaded_list = []
t_count = 0

# Open a logging file and begin the process of uploading
with open(where.LOG + '/Temp/' + temp_var['log_file_rank'], 'w') as log_w:
    for file in temp_var['small_list']:
        print('file: ', file)
        log = []
        log.append('file: %s' % file)
        # Create's or opens an existing server.  If the database cannot be created or opened it deletes and try again
        try:
            server = BioSeqDatabase.open_database(driver='sqlite3', db=ser_loc + '/' + temp_var['db_name'] + '.' + temp_var['key'] + '.db')
            print('server created')
            log.append('server created')
        except:
            print('server not created')
            log.append('server not created')
            os.remove(ser_loc + ('/%s.%s.db' % (temp_var['db_name'], temp_var['key'])))
            raise

        # Deprecated (all files are RNA, but I originally wanted to get the other types as well)
        s = str(file).lower()
        if s.find("rna") != -1:
            sub_db = 'RNA'
        elif s.find("protein") != -1:
            sub_db = 'PROTEIN'
        elif s.find("genomic") != -1:
            sub_db = 'GENOMIC'
        else:
            sub_db = 'Other'
        # Open the database while logging, and load the database with the file
        try:
            if sub_db not in server.keys():
                log.append('new database')
                print('new database')
                print('sub_db: ', sub_db)
                server.new_database(sub_db)
            db = server[sub_db]
            print('subdataase created')
            log.append('subdataase created')
            p = '/work5/r2294/DATA/PycharmProjects/GPCR_Orthologs/GPCR-Orthologs-Project/CODE/' \
                '1_Databases/NCBI_Data/refseq/release/multiprocessing/'

            count = db.load(SeqIO.parse(p + str(file), "genbank"))  # LOADING

            print('subdataase loaded')
            log.append('subdataase loaded')
            server.commit()
            print('Server Committed %s' % sub_db)
            print('%s database loaded with %s.' % (db.dbid, file))
            print("That file contains %s genbank records." % str(count))
            log.append('%s database loaded with %s.' % (db.dbid, file))
            t_count = t_count + count
            log_w.write('Process #%s.\n%s database loaded with %s, which contains %s GenBank files.  '
                        'The total number of files uploaded so far is %s.\n' % (process_number, db.dbid, file, count, t_count))
            log_w.write(str(log) + '\n')
            print('The total number of files loaded so far is %i.' % t_count)
            t = time.strftime("%I:%M:%S")
            loaded_list.append(("#%s#" % t) + ("<%s>  " % process_number) + str(file) + (' #GBK files: %s#' % count) + '\n')
        except:
            print('The database was unable to load from process %s' % process_number)
            log_w.write('The database was unable to load')
            log.append('The database was unable to load')
            log_w.write(str(log) + '\n')
            server.rollback()
            try:
                del server[sub_db]
                server.commit()
            except:
                raise
            raise
print("end of multi_dbupload.py from process %s" % process_number)
t = time.strftime("%I:%M:%S")
log_flag = True
while log_flag is True:
    try:
        with open(where.LOG + '/Ftp2Db_refseq-release-multiprocessing.log', 'a') as log_w:
            log_w.write('\n#%s# PROCESS %s uploaded %s GenBank files to %s.' % (t, process_number, t_count, (temp_var['db_name'] + '.' + temp_var['key'] + '.db')))
            for item in loaded_list:
                log_w.write(item)
        log_flag = False
    except:
        print('Could not open log_file from process %s' % process_number)
        time.sleep(30)
        log_flag = True
