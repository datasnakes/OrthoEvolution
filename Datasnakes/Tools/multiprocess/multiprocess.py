"""Use python's multiprocessing module to create multiple process and speed up
the completion of ClustalO and ETE3PAML.
"""
import re
import os
import fnmatch
import logging as log
import pandas as pd
from datetime import datetime as d
from multiprocessing import Pool
from time import time
import sys
from Datasnakes.Orthologs.Phylogenetics import ETE3PAML, RelaxPhylip
from Datasnakes.Orthologs.Align import ClustalO
from Datasnakes.Manager.utils import SplitList

## Home Directory
#home = os.getcwd() + '/'  # Use templating for directory name
#h = home
#os.chdir(h)
#
## Set up logging
#log.basicConfig(filename="kargpaml.log", level=log.INFO)
#log.info("#------------------------------------------------------------------")
#log.info("The script name is %s" % os.path.basename(__file__))
#log.info("The date and time is currently %s" % str(d.now()))
#log.info("#------------------------------------------------------------------")
#log.info("Run clustal omega.")
#
## Create the main output directories if they don't exist
#clustal_out = 'data/clustal-output/'
#paml_out = 'data/paml-output/'
#phyml_out = 'data/phyml-output/'
#
#dir_list = [clustal_out, phyml_out, paml_out]  # List of directories
#
#for directory in dir_list:
#    if os.path.exists(directory) is True:
#        log.info('The directory %s exists.' % directory)
#    else:
#        os.mkdir(directory)
#        log.info('The directory %s has been created' % directory)


def genes2analyze():
    """Get the list of genes based on pattern matching from the cds folder."""
    files = os.listdir('data/cds/')
    geneslist = []

    for filename in files:
        if fnmatch.fnmatch(filename, '*.ffn'):
            try:
                found = re.search('MASTER_(.+?)_CDS1.ffn', filename).group(1)
                geneslist.append(found)
            except AttributeError:
                log.error('There was an error in this attribute.')
        else:
            pass

    log.info("The genes are: %s " % geneslist)
    log.info("There are %s genes." % len(geneslist))

    # Save the list of genes that will be aligned & analyzed via PAML to the
    # clustal output dir
    df = pd.DataFrame(geneslist)
    df.to_csv(clustal_out + 'genes_to_align.txt', sep='\t',
              index=False, header=None)

    # Check to see if any alignment directories exist
    aligned = []
    for gene in geneslist:
        os.chdir(clustal_out)
        path = gene + "_Aligned/"
        if os.path.exists(path) is True:
            log.info("%s directory exists." % path)
            os.chdir(path)
            file = gene + "_aligned.phy"
            if os.path.exists(file) is True:
                log.info("%s exists." % file)
                aligned.append(gene)
                log.info("%s was removed from the list." % gene)
                os.chdir(h)
            else:
                os.chdir(h + clustal_out)
                os.system("rm -R " + path + " -f")
                log.info("The %s directory was deleted." % path)
                os.chdir(h)
        else:
            os.chdir(h)
            pass

    log.info("The aligned genes are: %s " % aligned)
    log.info("There are %s aligned genes." % len(aligned))

    # Create a new list of unaligned genes
    finalgeneslist = []
    for g in geneslist:
        if g not in aligned:
            finalgeneslist.append(g)

    log.info("The unaligned genes are: %s " % finalgeneslist)
    log.info("There are %s unaligned genes." % len(finalgeneslist))

    listgroups = SplitList(finalgeneslist, 'genes', n=int(30))
    aligneddict = SplitList(geneslist, 'genes', n=int(30))
    return listgroups, aligneddict


def multiclustal(gene):
    """Input a file of sequences to clustal omega in order to get a multiple
    sequnce alignment that will be analyzed for species divergence per gene
    with PAML."""
    # Create clustal omega gene directories
    gene_aligned_dir = clustal_out + gene + '_Aligned/'
    os.mkdir(gene_aligned_dir)

    # Copy file of sequences to output directory
    os.system('cp data/cds/MASTER_' + gene + '_CDS1.ffn ' + gene_aligned_dir)
    os.chdir(gene_aligned_dir)

    # Run clustal omega
    ClustalO(gene)
    log.info('Clustal omega has created a multiple sequence alignment for %s'
             % gene)

    # Run relaxphylip definition
    RelaxPhylip(gene)

    # Change to home directory
    os.chdir(h)


def multipaml(gene):
    """Input a phylip formated multiple sequence alignment to PAML in order to
    analyze divergence using the."""
    # Create paml gene directories
    gene_aligned_dir = clustal_out + gene + '_Aligned/'
    gene_paml_dir = paml_out + gene + '_PAML/'
    os.mkdir(gene_paml_dir)

    # Copy phylip file to directory
    os.system('cp ' + gene_aligned_dir + gene + '_aligned.phy ' +
              gene_paml_dir)
    os.chdir(gene_aligned_dir)

    ETE3PAML(gene)
    log.info('PAML has created output for %s' % gene)

    # Change to home directory
    os.chdir(h)


def main(function, geneslist):
    """This function uses a pool to start multiple processes to get clustal and
    PAML output. The argument (geneslist) should be a list of genes.
    """
    if len(geneslist) > 0:
        print(geneslist)
        ts = time()
        with Pool(processes=3) as p:
            p.map(function, geneslist)
            log.info("Multiprocessing with 8 processes has begun.")
            log.info("It took {} hours to get all algnments.".format((
                    time() - ts)/3600))
    elif len(geneslist) == 0:
        sys.exit("There are no genes that need aligned in your list.")
