Pandoc Documentation
====================

Use the ``docx2md.sh`` script to convert .docx files to .md (markdown)
format. The shell script uses pandoc to convert the files.

Dependencies
------------

`Pandoc <http://johnmacfarlane.net/pandoc/>`__ must be installed.

Setup
-----

**On Linux/Debian**

Make the script executable. Then run it. 1. ``chmod +x docx2md.sh`` 2.
``./docx2md.sh``

Examples
--------

In addition to a .sh/bash script to use with Pandoc, we've used
`pypandoc <>`__ to create a class that allows the conversion of
documents.

Convert markdown to docx
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools import PandocConverter
    PandocConverter(infile='README.md', outfmt='docx', outfile='README.docx')

Get a list of input formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools import PandocConverter
    PandocConverter.input_formats

    Out[17]:
    ['commonmark',
     'docbook',
     'docx',
     'epub',
     'haddock',
     'html',
     'json',
     'latex',
     'markdown',
     'markdown_github',
     'markdown_mmd',
     'markdown_phpextra',
     'markdown_strict',
     'mediawiki',
     'native',
     'odt',
     'opml',
     'org',
     'rst',
     't2t',
     'textile',
     'twiki']

Get a list of output formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    from OrthoEvol.Tools import PandocConverter
    PandocConverter.output_formats

    Out[18]:
    ['asciidoc',
     'beamer',
     'commonmark',
     'context',
     'docbook',
     'docbook5',
     'docx',
     'dokuwiki',
     'dzslides',
     'epub',
     'epub3',
     'fb2',
     'haddock',
     'html',
     'html5',
     'icml',
     'json',
     'latex',
     'man',
     'markdown',
     'markdown_github',
     'markdown_mmd',
     'markdown_phpextra',
     'markdown_strict',
     'mediawiki',
     'native',
     'odt',
     'opendocument',
     'opml',
     'org',
     'plain',
     'revealjs',
     'rst',
     'rtf',
     's5',
     'slideous',
     'slidy',
     'tei',
     'texinfo',
     'textile',
     'zimwiki']

