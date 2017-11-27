"""Class to access NCBI's ftp server and easily download databases."""
import tarfile
import re
from time import time
from datetime import datetime
from multiprocessing.pool import ThreadPool
import os
from shutil import make_archive
from ftplib import error_perm, all_errors
# from progress.bar import Bar
# TODO Create a progress bar; Integrate with Threading/downloading

from OrthoEvol.Tools.ftp.baseftp import BaseFTPClient
from OrthoEvol.Tools.logit import LogIt


class NcbiFTPClient(BaseFTPClient):
    """Access NCBI's FTP servers with ease."""

    def __init__(self, email):
        _ncbi = 'ftp.ncbi.nlm.nih.gov'
        super().__init__(_ncbi, email, keepalive=False, debug_lvl=0)
        self._datafmt = '%m-%d-%Y@%I:%M:%S-%p'
        self._date = str(datetime.now().strftime(self._datafmt))
        self.blastpath = '/blast/'
        self.blastdb_path = '/blast/db/'
        self.blastfasta_path = '/blast/db/FASTA/'
        self.refseqrelease_path = '/refseq/release/'
        self.windowmasker_path = self.blastpath + 'windowmasker_files/'
        self._taxdb = 'taxdb.tar.gz'  # Located in self.blastdb_path

        # TODO Use Turn into a json file, dict, or config
        self.blastdbs = []
        self.blastfastadbs = []

        # TODO Create dictionary of refseqrelease dbs, seqtypes, filetypes
        self.refseqreleasedbs = []
        self.refseqrelease_seqtypes = []
        self.refseqrelease_filetypes = []

        # Set up logger
        self.ncbiftp_log = LogIt().default(logname="NCBI-FTP", logfile=None)

    @classmethod
    def _pathformat(cls, path):
        """Ensure proper formatting of the path.

        :param path: FTP path.
        :type path: str
        """
        pattern = re.compile('^/(.*?)/$')
        if not re.match(pattern, path):
            raise ValueError('Your path is not in a proper format.')

    @classmethod
    def _archive(cls, archive_name, folder2archive, archive_type='gztar'):
        """Archive all the files in the folder and compress the archive.

        :param archive_name:
        :param folder2archive:
        :param archive_type:
        :return:
        """
        os.chdir(folder2archive)  # Enter the full path
        os.chdir('..')
        archive_location = os.path.join(os.getcwd(), archive_name)
        os.chdir(folder2archive)
        make_archive(archive_location, folder2archive, archive_type)

    def walk(self, path):
        """Walk the ftp server and get files and directories.

        :param path: FTP path to be walked.
        """
        file_list, dirs, nondirs = [], [], []
        try:
            self.ftp.cwd(path)
        except error_perm as ep:
            self.ncbiftp_log.info("Current path: %s" % self.ftp.pwd() + ep.__str__() + path)
            return [], []
        else:
            self.ftp.retrlines('LIST', lambda x: file_list.append(x.split()))
            for info in file_list:
                ls_type, name = info[0], info[-1]
                if ls_type.startswith('d'):
                    dirs.append(name)
                else:
                    nondirs.append(name)
            return dirs, nondirs

    def download_file(self, filename):
        """Download the files one by one.

        :param filename:
        :return:
        """
        with open(filename, 'wb') as localfile:
            self.ftp.retrbinary('RETR %s' % filename, localfile.write)
            self.ncbiftp_log.info('%s was downloaded.' % str(filename))

    def _download_windowmasker(self, windowmaskerfilepath):
        """Download the window masker files.

        :param windowmaskerfilepath: FTP path to windowmasker files.
        """
        wmsplit = windowmaskerfilepath.split(sep='/')
        wmdir = wmsplit[0]
        wmfile = wmsplit[1]

        os.makedirs(wmdir, exist_ok=True)

        if not os.path.exists(os.path.join(wmdir, wmfile)):
            try:
                with open(os.path.join(wmdir, wmfile), 'wb') as localfile:
                    self.ftp.retrbinary('RETR %s' % windowmaskerfilepath, localfile.write)
                    self.ncbiftp_log.info('%s was downloaded.' % str(windowmaskerfilepath))
            except all_errors:
                os.remove(os.path.join(wmdir, wmfile))
                err_msg = '%s could not be download and was deleted.'
                self.ncbiftp_log.error(err_msg % str(windowmaskerfilepath))
        else:
            self.ncbiftp_log.info('%s exists.' % str(windowmaskerfilepath))

    @classmethod
    def extract_file(cls, file2extract):
        """Extract a tar.gz file.

        :param file2extract: Path to the tar.gz file to extract.
        """
        if str(file2extract).endswith('tar.gz'):
            tar = tarfile.open(file2extract)
            tar.extractall()
            tar.close()
            log_msg = 'Files were successfully extracted from %s'
            cls.ncbiftp_log.info(log_msg % file2extract)

    def listfiles(self, path='cwd'):
        """List all files in a path.

        :param path: Directory path.
        :return: A list of files in a directory.
        """
        if path == 'cwd':
            path = self.ftp.pwd()
        path = path
        self._pathformat(path)
        _, files = self.walk(path)
        return files

    def listdirectories(self, path='cwd'):
        """List all directories in a path.

        :param path: Directory path.
        :return: A list of subdirectories.
        """
        if path == 'cwd':
            path = self.ftp.pwd()
        path = path
        self._pathformat(path)
        directories, _ = self.walk(path)
        return directories

    def getwindowmaskerfiles(self, taxonomy_ids, download_path):
        """Download NCBI's window masker binary files for each taxonomy id.

        :param taxonomy_ids:
        :param download_path:
        :return:
        """
        self.ftp.cwd(self.windowmasker_path)
        taxonomy_dirs = self.listdirectories(self.windowmasker_path)

        # Change to directory input
        unfound_species = []
        found_species = []
        for taxonomy_id in taxonomy_ids:
            if str(taxonomy_id) not in taxonomy_dirs:
                self.ncbiftp_log.warning('%s does not exist.' % taxonomy_id)
                unfound_species.append(taxonomy_id)
            else:
                found_species.append(taxonomy_id)

        windowmaskerfiles = []
        for found_taxid in found_species:
            filepath = str(found_taxid) + '/wmasker.obinary'
            windowmaskerfiles.append(filepath)

        self.ncbiftp_log.info('You are about to download theses files: %s' %
                              windowmaskerfiles)

        # Move to directory for file downloads
        os.chdir(download_path)

        # Download the files using multiprocessing
        download_time_secs = time()
        with ThreadPool(1) as download_pool:
            download_pool.map(self._download_windowmasker, windowmaskerfiles)
            minutes = round(((time() - download_time_secs) / 60), 2)
        self.ncbiftp_log.info("Took %s minutes to download the files." %
                              minutes)
        return unfound_species

    def getrefseqrelease(self, taxon_group, seqtype, seqformat, download_path,
                         extract=True):
        """Download the refseq release database.

        :param taxon_group:
        :param seqtype:
        :param seqformat:
        :param download_path:
        :param extract:
        :return:
        """
        self.ftp.cwd(self.refseqrelease_path)
        taxon_dirs = self.listdirectories(self.refseqrelease_path)

        # Change to directory input
        if taxon_group not in taxon_dirs:
            raise FileNotFoundError('%s does not exist.' % taxon_group)

        self.ftp.cwd(taxon_group)
        curpath = self.ftp.pwd() + '/'
        releasefiles = self.listfiles(curpath)

        files2download = []
        pattern = re.compile('^' + taxon_group + '[.](.*?)[.]' + seqtype
                             + '[.]' + seqformat + '[.]gz$')
        for releasefile in releasefiles:
            if re.match(pattern, releasefile):
                files2download.append(releasefile)

        self.ncbiftp_log.info('You are about to download theses files: %s' %
                              files2download)

        # Move to directory for file downloads
        os.chdir(download_path)

        # Download the files using multiprocessing
        download_time_secs = time()
        with ThreadPool(1) as download_pool:
            download_pool.map(self.download_file, files2download)
            minutes = round(((time() - download_time_secs) / 60), 2)
        self.ncbiftp_log.info("Took %s minutes to download the files." %
                              minutes)

        if extract:
            extract_time_secs = time()
            with ThreadPool(1) as extract_pool:
                extract_pool.map(self.extract_file, files2download)
                minutes = round(((time() - extract_time_secs) / 60), 2)
            self.ncbiftp_log.info("Took %s minutes to extract from all files." %
                                  minutes)

    def getblastfasta(self, database_name, download_path, extract=True):
        """Download the fasta sequence database (not formatted).

        :param database_name:
        :param download_path:
        :param extract:
        """
        if str(database_name).startswith('est'):
            raise NotImplementedError('Est dbs cannot be downloaded yet.')
        self.ftp.cwd(self.blastfasta_path)
        blastfastafiles = self.listfiles(self.blastfasta_path)

        files2download = []
        for dbfile in blastfastafiles:
            if database_name in str(dbfile):
                files2download.append(dbfile)

        # Append the taxonomy database
        files2download.append(self._taxdb)
        self.ncbiftp_log.info('You are about to download theses files: %s' %
                              files2download)

        # Move to directory for file downloads
        os.chdir(download_path)

        # Download the files using multiprocessing
        download_time_secs = time()
        with ThreadPool(1) as download_pool:
            download_pool.map(self.download_file, files2download)
            minutes = round(((time() - download_time_secs) / 60), 2)
        self.ncbiftp_log.info("Took %s minutes to download the files." %
                              minutes)

        if extract:
            extract_time_secs = time()
            with ThreadPool(1) as extract_pool:
                extract_pool.map(self.extract_file, files2download)
                minutes = round(((time() - extract_time_secs) / 60), 2)
            self.ncbiftp_log.info("Took %s minutes to extract from all files." %
                                  minutes)

    def getblastdb(self, database_name, download_path, extract=True):
        """Download the formatted blast database.

        :param database_name: Name of preformatted blastdb to download.
        :param download_path: Directory path/location to download blastdb.
        :param extract (bool): True or False for extract tar.gz db files.
        """
        if str(database_name).startswith('est'):
            raise NotImplementedError('Est dbs cannot be downloaded yet.')
        self.ftp.cwd(self.blastdb_path)
        blastdbfiles = self.listfiles(self.blastdb_path)

        files2download = []
        for dbfile in blastdbfiles:
            if database_name in str(dbfile):
                files2download.append(dbfile)

        # Append the taxonomy database
        files2download.append(self._taxdb)

        # Move to directory for file downloads
        os.chdir(download_path)

        absentfiles = []

        # Ensure that files aren't already downloaded
        for file2download in files2download:
            if not os.path.exists(os.path.join(download_path, file2download)):
                absentfiles.append(file2download)

        if len(absentfiles) > 0:
            self.ncbiftp_log.info('You are about to download these files: %s\n' %
                                  absentfiles)
            # Download the files using multiprocessing
            download_time_secs = time()
            with ThreadPool(1) as download_pool:
                download_pool.map(self.download_file, files2download)
                minutes = round(((time() - download_time_secs) / 60), 2)
            self.ncbiftp_log.info("Took %s minutes to download the files.\n" %
                                  minutes)

        if extract:
            self.ncbiftp_log.info('Now it\'s time to extract files.')
            extract_time_secs = time()
            with ThreadPool(3) as extract_pool:
                extract_pool.map(self.extract_file, files2download)
                minutes = round(((time() - extract_time_secs) / 60), 2)
            self.ncbiftp_log.info("Took %s minutes to extract from all files.\n" %
                                  minutes)

        # Remove all tar.gz files
        curfiles = os.listdir()
        for curfile in curfiles:
            if str(curfile).endswith('tar.gz'):
                os.remove(curfile)

    def updatedb(self, database_path=os.getcwd(), update_days=7):
        """Check for when the database was last updated.

        Refseq release databases should only be updated every few months.

        :param database_path: Directory path of existing database.
        :param update_days (int): Number of days to update.
        :return:
        """
        # TODO Prevent users from updated refseq if certain days
        # Get a list of the files in the path
        filesinpath = os.listdir(database_path)
        for fileinpath in filesinpath:
            if str(fileinpath).endswith('.nal'):
                nalfile = str(fileinpath)
                dbname, ext = nalfile.split('.')
                filetime = datetime.fromtimestamp(os.path.getctime(nalfile))
                format_filetime = filetime.strftime("%b %d, %Y at %I:%M:%S %p")

            elif str(fileinpath).endswith('.gbff'):
                gbff_file = str(fileinpath)
                taxon_group, _, seqtype, ext = gbff_file.split('.')
                filetime = datetime.fromtimestamp(os.path.getctime(gbff_file))
                format_filetime = filetime.strftime("%b %d, %Y at %I:%M:%S %p")

        self.ncbiftp_log.info("Your database was last updated on: %s" %
                              format_filetime)

        time_elapsed = datetime.now() - filetime

        if ext == 'nal' and time_elapsed.days >= update_days:
            self.ncbiftp_log.warning('\nYour blast database needs updating.')
            archive_name = "blastdb_archive_" + self._date
            self._archive(archive_name, folder2archive=database_path)
            self.getblastdb(dbname, download_path=database_path, extract=True)

        elif ext == 'gbff' and time_elapsed.days >= 70:
            self.ncbiftp_log.warning('\nYour refseq release database needs updating.')
            archive_name = "refseqrelease_archive_" + self._date
            self._archive(archive_name, folder2archive=database_path)

            # TODO Create a way to get seqtype and seqformat from filename
            self.getrefseqrelease(taxon_group, seqtype, ext, database_path,
                                  extract=True)

        # TODO Add elif for handling databases not handled by this class.

        elif ext != 'nal' or 'gbff':
            raise NotImplementedError("Updating for that database is unsupported.")

        else:
            self.ncbiftp_log.info('\nYour database is still up to date.')
