"""Use python's multiprocessing module to create multiple process and speed up
the completion of classes."""
import re
import os
import fnmatch
import pandas as pd
from multiprocessing import Pool
from time import time
import sys
from Datasnakes.Orthologs.Phylogenetics import ETE3PAML, RelaxPhylip
from Datasnakes.Orthologs.Align import ClustalO
from Datasnakes.Manager.utils import SplitList


class MultiProc(object):
    def __init__(self, workingdir=os.getcwd()):
        """Set up Align and phylogenetics data directories."""
        self.h = workingdir
        os.chdir(self.h)
        # Create the main output directories if they don't exist
        self.clustal_out = 'data/clustal-output/'
        self.paml_out = 'data/paml-output/'
        self.phyml_out = 'data/phyml-output/'

        # List of directories
        dir_list = [self.clustal_out, self.phyml_out, self.paml_out]

        for directory in dir_list:
            if os.path.exists(directory) is True:
                print('The directory %s exists.' % directory)
            else:
                os.mkdir(directory)
                print('The directory %s has been created' % directory)

    def genes2analyze(self):
        """Get a list of genes based on files in the cds folder."""
        files = os.listdir('data/cds/')
        geneslist = []

        for filename in files:
            if fnmatch.fnmatch(filename, '*.ffn'):
                try:
                    found = re.search('MASTER_(.+?)_CDS1.ffn',
                                      filename).group(1)
                    geneslist.append(found)
                except AttributeError:
                    print('There was an error in this attribute.')
            else:
                pass

        print("The genes are: %s " % geneslist)
        print("There are %s genes." % len(geneslist))

        # Save the list of genes that will be aligned & analyzed via PAML to
        # the clustal output dir
        df = pd.DataFrame(geneslist)
        df.to_csv(self.clustal_out + 'genes_to_align.txt', sep='\t',
                  index=False, header=None)

        # Check to see if any alignment directories exist
        aligned = []
        for gene in geneslist:
            os.chdir(self.clustal_out)
            path = gene + "_Aligned/"
            if os.path.exists(path) is True:
                print("%s directory exists." % path)
                os.chdir(path)
                file = gene + "_aligned.phy"
                if os.path.exists(file) is True:
                    print("%s exists." % file)
                    aligned.append(gene)
                    print("%s was removed from the list." % gene)
                    os.chdir(self.h)
                else:
                    os.chdir(os.path.join(self.h, self.clustal_out))
                    os.system("rm -R " + path + " -f")
                    print("The %s directory was deleted." % path)
                    os.chdir(self.h)
            else:
                os.chdir(self.h)
                pass

        print("The aligned genes are: %s " % aligned)
        print("There are %s aligned genes." % len(aligned))

        # Create a new list of unaligned genes
        finalgeneslist = []
        for g in geneslist:
            if g not in aligned:
                finalgeneslist.append(g)

        print("The unaligned genes are: %s " % finalgeneslist)
        print("There are %s unaligned genes." % len(finalgeneslist))

        listgroups = SplitList(finalgeneslist, 'genes', n=int(30))
        aligneddict = SplitList(geneslist, 'genes', n=int(30))
        return listgroups, aligneddict

    def multiclustal(self, gene):
        """Multiprocessing function with ClustalO"""
        # Create clustal omega gene directories
        gene_aligned_dir = self.clustal_out + gene + '_Aligned/'
        os.mkdir(gene_aligned_dir)

        # Copy file of sequences to output directory
        os.system('cp data/cds/MASTER_' + gene + '_CDS1.ffn ' +
                  gene_aligned_dir)
        os.chdir(gene_aligned_dir)

        # Run clustal omega
        ClustalO(gene)
        print('Clustal omega has created a multiple sequence alignment for %s'
              % gene)

        # Run relaxphylip definition
        RelaxPhylip(gene)

        # Change to home directory
        os.chdir(self.h)

    def multipaml(self, gene):
        """Input a phylip formated multiple sequence alignment to PAML in order to
        analyze divergence using the."""
        # Create paml gene directories
        gene_aligned_dir = self.clustal_out + gene + '_Aligned/'
        gene_paml_dir = self.paml_out + gene + '_PAML/'
        os.mkdir(gene_paml_dir)

        # Copy phylip file to directory
        os.system('cp ' + gene_aligned_dir + gene + '_aligned.phy ' +
                  gene_paml_dir)
        os.chdir(gene_aligned_dir)

        ETE3PAML(gene)
        print('PAML has created output for %s' % gene)

        # Change to home directory
        os.chdir(self.h)

    def main(function, geneslist, processes):
        """This function uses a pool to start multiple processes to get clustal and
        PAML output. The argument (geneslist) should be a list of genes.
        """
        if len(geneslist) > 0:
            print(geneslist)
            ts = time()
            with Pool(processes) as p:
                p.map(function, geneslist)
                print("Multiprocessing with 8 processes has begun.")
                print("It took {} hours to get all algnments.".format((
                        time() - ts)/3600))
        elif len(geneslist) == 0:
            sys.exit("There are no genes that need aligned in your list.")
