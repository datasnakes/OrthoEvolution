from shutil import copyfile, move
import os
from os.path import exists, join
import pypandoc

from Datasnakes.Tools.logit import LogIt


class PandocConverter(object):
    """Convert files using the pypandoc wrapper.

    :param infile:  path to the input file
    """
    def __init__(self, infile, outfile_path):
        self.infile = infile
        self.outfile_path = outfile_path
        self.pandoc_log = LogIt().default('pandoc converter', None)
        self.md2rst()

    def md2rst(self):
        infile = str(self.infile)
        filepath, ext = infile.split('.')
        splitfilepath = filepath.split(os.sep)
        filename = splitfilepath[-1]
        outfile = join(self.outfile_path, filename + ".rst")
        if ext == 'md':
            print(outfile)
            pypandoc.convert_file(infile, 'rst',  outputfile=outfile)
            self.pandoc_log.info('%s.md converted to %s.rst.' % (filename, filename))
        else:
            self.pandoc_log.error('%s not supported by this function.' % ext)

    def archivemdfiles(self):
        if not exists('archive'):
            os.mkdir('archive')
            self.pandoc_log.info("The archive directory was created.")

        move(self.infile, join("archive", self.infile))
        self.pandoc_log.info("%s was moved to archive directory." % self.infile)


class CreateDocs(object):
    def __init__(self):
        self.createdocs_log = LogIt().default('create docs', None)
        self.current_dir = os.path.dirname(os.path.realpath(__file__))
        self._docsdir = os.path.join(self.current_dir, 'docs')
        self.up = '..'
        self.docs_source = os.path.join(self._docsdir, 'source')
        self.main_readme = self._pathtomainreadme()
        self.readme2index()
        self.packagename = 'Datasnakes'

    def readme2index(self):
        """Copy README.rst from the top level of your repo to docs source."""
        indexpath = os.path.join(self.docs_source, 'index.rst')
        copyfile(self.main_readme, indexpath)

    def _pathtomainreadme(self):
        os.chdir(self.current_dir)
        os.chdir(self.up)
        readmedir = os.getcwd()
        main_readme = os.path.join(readmedir, 'README.rst')
        os.chdir(self.current_dir)
        return main_readme

    def getfiles2convert(self):
        """Get a list of files to convert."""
        os.chdir(self.current_dir)
        os.chdir(self.up)
        datasnakesdir = os.path.join(os.getcwd(), self.packagename)

        skip = ['new_basic_project', 'new_repository', 'new_database_repo',
                'index', 'new_research', 'new_user', 'web', 'data',
                'new_website', 'config']

        files2convert = []
        for dirname, subdirlist, filelist in os.walk(datasnakesdir):
            subdirlist[:] = [d for d in subdirlist if d not in skip]
            for filename in filelist:
                if filename == 'README.md':
                    file = os.path.join(dirname, filename)
                    files2convert.append(file)

        return files2convert


if __name__ == "__main__":
    createdocs = CreateDocs()
    files2convert = createdocs.getfiles2convert()
    createdocs.readme2index()

    for file2convert in files2convert:
        PandocConverter(file2convert, outfile_path=createdocs.docs_source)
