#!/usr/bin/env python
"""This standalone script downloads files from NCBI's ftp."""
import os
import fnmatch
from subprocess import call, CalledProcessError
import contextlib
import argparse
import textwrap
import sys

from OrthoEvol.Tools.ftp import NcbiFTPClient

# Raise an error if you're not on linux. Windows generally doesn't have wget.
if sys.platform != 'linux':
    msg = 'This interface is not intended for use on your platform.'
    raise NotImplementedError(msg)


def write_to_file(hostname, dbname, dbpath, filenames):
    # Create a for loop that writes the list/text file of files wanted
    with open('downloadlist.txt', 'w') as downloads:
        for filename in filenames:
            # Get only those files.
            if fnmatch.fnmatch(filename, dbname + '*'):
                refseq_file = os.path.join(filename)
                # Write the url of each refseq_rna db file to a text file.
                downloads.writelines(hostname + dbpath + refseq_file + '\n')
            # use elif here to get the taxdb.tar.gz file.
            elif fnmatch.fnmatch(filename, 'taxdb*'):
                taxdb_file = os.path.join(filename)
                downloads.writelines(hostname + dbpath + taxdb_file + '\n')


def main(email, dbtype, dbname, num_procs=8):
    """[summary]

    :param email: [description]
    :type email: [type]
    :param dbtype: [description]
    :type dbtype: [type]
    :param dbname: [description]
    :type dbname: [type]
    :param num_procs: The number of processors to use, defaults to 8
    :type num_procs: int, optional
    :raises NotImplementedError: [description]
    """
    ncbiftp = NcbiFTPClient(email=email)
    log = ncbiftp.ncbiftp_log

    if dbtype == 'blastdbv5':
        # This is a list of the file names in the current directory
        dbpath = ncbiftp.blastdbv5_path
        filenames = ncbiftp.listfiles(dbpath)

        write_to_file(ncbiftp.ftp.host, dbname, dbpath, filenames)

    elif dbtype == 'blastdb':
        # This is a list of the file names in the current directory
        dbpath = ncbiftp.blastdb_path
        filenames = ncbiftp.listfiles(dbpath)

        write_to_file(ncbiftp.ftp.host, dbname, dbpath, filenames)
    else:
        raise NotImplementedError('That database is unsupported.')

    # Download the list of files using 'wget' on linux/unix
    with contextlib.suppress(os.error):
        cmd = 'cat downloadlist.txt | xargs -n 1 -P ' + \
            int(num_procs) + ' wget'
        status = call([cmd], shell=True)
        if status == 0:
            log.info("The %s %s files have downloaded." % (dbname, dbtype))
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
    parser.add_argument('-dbname', '--database-name', help='The name or seqtype of the database',
                        required=True)
    parser.add_argument('-n', '--num-procs',
                        help='Enter the number of processors to use to download the files',
                        required=False)

    args = parser.parse_args()

    main(args.email, args.dbtype, args.dbname, args.num_procs)
