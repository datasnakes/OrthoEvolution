import os
from os.path import isfile, join, basename, splitext, exists
import psutil
from subprocess import PIPE
from shutil import copyfile


def md2rst(path):
    """Converts all markdown documents to rst documents in a folder.

    It does not include the readme, but it does copy the readme and rename it
    as the index to be used in the docs for a folder.
    """
    home = path  # Use this directory
    os.chdir(home)
    readme = 'README.md'
    copyfile(readme, 'index.md')  # Copy the readme as an index file.
    exclude = [readme]  # Exclude the readme
    files = [f for f in os.listdir(home) if isfile(join(home, f))]

    # Use a loop to convert the docs.
    for file in files:
        if file not in exclude:  # Include all files except 'README.md'
            if file.endswith('.md'):  # Only convert markdown files
                name = splitext(basename(file))[0]  # Get the file name w/o ext
                cmd = ["pandoc", "--from=markdown", "--to=rst", "--output=" +
                       name + ".rst", name + ".md"]
                convert = psutil.Popen(cmd, stdout=PIPE)
                convert.wait()
                print("%s was converted to rst format." % file)

                # If an archive directory doesn't exist, create
                # If it does exist, markdown files will be moved there.
                if not exists('archive'):
                    os.mkdir('archive')
                    print("The archive directory was created.")
                cmd2 = ["mv", name + ".md", "archive/"]
                move = psutil.Popen(cmd2, stdout=PIPE)
                move.wait()
                print("%s was moved to archive directory." % file)
