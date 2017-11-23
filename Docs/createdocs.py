#from shutil import copyfile
#import os
#from os.path import join
#import configparser
#
#config = configparser.ConfigParser()
#
#config.read('docs.cfg')
#
#try:
#    import pypandoc
#except Exception:
#    from pypandoc.pandoc_download import download_pandoc
#    download_pandoc()
#
#from OrthoEvol.Tools.logit import LogIt
#
#
#class PandocConverter(object):
#    """Convert files using the pypandoc wrapper.
#
#    :param infile:  path to the input file
#    :param outfile_path:  output path
#    """
#    def __init__(self, infile, outfile_path):
#        self.infile = infile
#        self.outfile_path = outfile_path
#        self.pandoc_log = LogIt().default('pandoc converter', None)
#        self.md2rst()
#
#    def md2rst(self):
#        infile = str(self.infile)
#        filepath, ext = infile.split('.')
#        splitfilepath = filepath.split(os.sep)
#        initial_filename = splitfilepath[-1]
#        final_filename = str(splitfilepath[-2] + initial_filename).lower()
#        outfile = join(self.outfile_path, final_filename + ".rst")
#        if ext == 'md':
#            pypandoc.convert_file(infile, 'rst',  outputfile=outfile)
#            self.pandoc_log.info('%s.md converted to %s.rst.' %
#                                 (infile, final_filename))
#        else:
#            self.pandoc_log.error('%s not supported by this function.' % ext)
#
#
#class CreateDocs(object):
#    def __init__(self):
#        """"Initialize this class."""
#        self.createdocs_log = LogIt().default('create docs', None)
#        self.current_dir = os.path.dirname(os.path.realpath(__file__))
#        self._docsdir = os.path.join(self.current_dir, 'docs')
#        self.up = '..'
#        self.docs_source = os.path.join(self._docsdir, 'source')
#        self.packagename = 'OrthoEvol'
#        self.convertfiles()
#
#    def docscfg2lists(self):
#        config = configparser.ConfigParser()
#        os.chdir(self.current_dir)
#        config.read('docs.cfg')
#        sections = config.sections()
#        cookiesfiles = []
#        managerfiles = []
#        orthologsfiles = []
#        pipelinefiles = []
#        toolsfiles = []
#        for section in sections:
#            for key in config[section]:
#                if section == 'cookies':
#                    cookiesfiles.append(config[section][key])
#                elif section == 'manager':
#                    managerfiles.append(config[section][key])
#                elif section == 'orthologs':
#                    orthologsfiles.append(config[section][key])
#                elif section == 'pipeline':
#                    pipelinefiles.append(config[section][key])
#                elif section == 'tools':
#                    toolsfiles.append(config[section][key])
#        return (cookiesfiles, managerfiles, orthologsfiles, pipelinefiles,
#                toolsfiles)
#
#    def getfiles2convert(self):
#        """Get a list of files to convert."""
#        os.chdir(self.current_dir)
#        os.chdir(self.up)
#        datasnakesdir = os.path.join(os.getcwd(), self.packagename)
#
#        # Skip these folders when looking for README.md
#        # Most of these are folders in the cookiecutter repository
#        skip = ['new_basic_project', 'new_repository', 'new_database_repo',
#                'index', 'new_research', 'new_user', 'web', 'data',
#                'new_website', 'config', 'utils']
#
#        files2convert = []
#
#        # Use os.walk to walk from the datasnakes/package directory and down
#        # to find README.md files in each main submodule.
#        for dirname, subdirlist, filelist in os.walk(datasnakesdir):
#            subdirlist[:] = [d for d in subdirlist if d not in skip]
#            for filename in filelist:
#                if filename == 'README.md':
#                    file = os.path.join(dirname, filename)
#                    files2convert.append(file)
#
#        return files2convert
#
#    def convertfiles(self):
#        """Convert a list of files from .md to .rst."""
#        files2convert = self.getfiles2convert()
#        cookiesfiles, managerfiles, orthologsfiles, pipelinefiles, toolsfiles = self.docscfg2lists()
#        for file2convert in files2convert:
#            splitfilepath = str(file2convert).split(os.sep)
#            filename = str(splitfilepath[-2]).lower() + 'readme.rst'
#            if filename in cookiesfiles:
#                outpath = join(self.docs_source, 'cookies')
#                PandocConverter(file2convert, outpath)
#            elif filename in managerfiles:
#                outpath = join(self.docs_source, 'manager')
#                PandocConverter(file2convert, outpath)
#            elif filename in orthologsfiles:
#                outpath = join(self.docs_source, 'orthologs')
#                PandocConverter(file2convert, outpath)
#            elif filename in pipelinefiles:
#                outpath = join(self.docs_source, 'pipeline')
#                PandocConverter(file2convert, outpath)
#            elif filename in toolsfiles:
#                outpath = join(self.docs_source, 'tools')
#                PandocConverter(file2convert, outpath)
#            elif filename == 'orthoevolreadme.rst':
#                outpath = join(self.docs_source, 'tutorial')
#                PandocConverter(file2convert, outpath)
#            else:
#                pass
#        self.createdocs_log.info('Your docs were converted to .rst format.')
#
#
#if __name__ == "__main__":
#    CreateDocs()
