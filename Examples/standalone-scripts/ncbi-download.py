"""This standalone script downloads files from NCBI's ftp."""
from OrthoEvol.Tools.ftp import NcbiFTPClient
import os
import fnmatch
from subprocess import call, CalledProcessError
import contextlib
import argparse
import textwrap
import sys

# Raise an error if you're not on linux. Windows generally doesn't have wget.
if 'linux' not in str(sys.platform):
    msg = 'This interface is not intended for use on your platform.'
    raise NotImplementedError(msg)


def main(email, dbtype, dbname, preformatted, num_procs=8):
    ncbiftp = NcbiFTPClient(email=email)
    log = ncbiftp.ncbiftp_log
    accepted = ['yes', 'Yes', 'y', 'Y']


    if dbtype == 'blastdb' and preformatted in accepted:

        # This is a list of the file names in the current directory
        dbpath = ncbiftp.blastdb_path
        filenames = ncbiftp.listfiles(dbpath)

        # Create a for loop that writes the list/text file of files wanted
        with open('downloadlist.txt', 'w') as downloads:
            for filename in filenames:
                if fnmatch.fnmatch(filename, dbname + '*'):  # Get only those files.
                    refseq_file = os.path.join(filename)
                    # Write the url of each refseq_rna db file to a text file.
                    downloads.writelines(ncbiftp.ftp.host + dbpath + refseq_file + '\n')
                # use elif here to get the taxdb.tar.gz file.
                elif fnmatch.fnmatch(filename, 'taxdb*'):
                    taxdb_file = os.path.join(filename)
                    downloads.writelines(ncbiftp.ftp.host + dbpath + taxdb_file + '\n')

    elif preformatted not in accepted:
        raise NotImplementedError('Non-formatted databases are NOT unsupported.')

    else:
        raise NotImplementedError('That database is unsupported.')

    # Download the list of files using 'wget' on linux/unix
    with contextlib.suppress(os.error):
        cmd = 'cat downloadlist.txt | xargs -n 1 -P ' + int(num_procs) + ' wget'
        status = call([cmd], shell=True)
        if status == 0:
            log.info("The %s blast db files have downloaded." % dbname)
        else:
            log.error(CalledProcessError)
            ncbiftp.close_connection()

    ncbiftp.close_connection()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                    This is a command line interface to download
                                    NCBI Blast Databases.

                                    Preformatted blast databases are currently the only
                                    supported databases '''))
    parser.add_argument('-e', '--email', help='Enter your email address',
                        required=True)
    parser.add_argument('-dbtype', '--database-type',
                        help='Enter the name of the NCBI database.',
                        required=True)
    parser.add_argument('-dbname', '--database-name', help='Respond with yes or no',
                        required=True)
    parser.add_argument('-p', '--preformatted', help='Respond with yes or no',
                        required=True)
    parser.add_argument('-n', '--num-procs',
                        help='Enter the number of processors to use to download the files',
                        required=False)

    args = parser.parse_args()

    main(args.email, args.dbtype, args.dbtype, args.preformatted, args.num_procs)