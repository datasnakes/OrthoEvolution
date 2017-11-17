# -*- coding: utf-8 -*-
import pypandoc

from Datasnakes.Tools import LogIt


class PandocConverter(object):
    input_formats, output_formats = pypandoc.get_pandoc_formats()
    """Convert files using the pypandoc wrapper.

    :param infile:  path to the input file
    :param outfile_path:  output path
    """
    def __init__(self, infile, outfmt, outfile):
        self.infile = infile
        self.outfmt = outfmt
        self.outfile = outfile
        self.pandoc_inputs, self.pandoc_outputs = pypandoc.get_pandoc_formats()
        self.pandoc_log = LogIt().default('pandoc converter', None)
        self.convert()

    def convert(self):
        infile = str(self.infile)
        filename, ext = infile.split('.')
        try:
            pypandoc.convert_file(infile, self.outfmt, self.outfile)
            self.pandoc_log.info('%s converted to %s' %
                                 (infile, self.outfile))
        except Exception as err:
            self.pandoc_log.error(err)
