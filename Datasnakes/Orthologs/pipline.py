from Datasnakes.Orthologs.Align.QualityControl.filter import FilteredAlignment
from Datasnakes.Orthologs.Phylogenetics.IQTree.best_tree import FilteredTree
from Datasnakes.Orthologs.Phylogenetics.PAML.codeml import CodemlRun
import os
import time
from shutil import copy
from pathlib import Path


class OrthologPipeline(object):

    def __init__(self, genes, qsub_template, python_file, home=os.getcwd()):
        self.index = Path(home)  # index directory
        self.raw_data = self.index.parent / Path('raw_data')
        self.data = self.index.parent / Path('data')
        with open(genes, 'r') as gene_list:
            for gene in gene_list.readlines():
                qsub_file = copy(str(self.index / Path(qsub_template)),
                                 str(self.raw_data / Path(gene) / Path(gene + '_job.sh')))

                self.na_fasta = str(self.raw_data / Path(gene + '.ffn'))
                self.aa_fasta = str(self.raw_data / Path(gene + '.faa'))
                self.P2N_na_alignment = str(self.raw_data / Path(gene + '_P2N_na.aln'))
                self.iqtree_newick = str(self.raw_data / Path(gene + 'iqtree.nwk'))
                self.paml_control = str(self.raw_data / Path('paml.ctl'))


    def batch_script_setup(self, qsub_file, python_file, raw_data_path, author, email, memory, gene):
        date_time = time.strftime("%H:%M:%S on %d/%m/%Y")
        qsub_command = "qsub -v PYTHONFILE=%s RAWDATA=%s AUTHOR=%s EMAIL=%s GBS=%s GENE=%s DATE=%s %s" % \
                       (python_file, raw_data_path, author, email, memory, gene, date_time, qsub_file)
    # qca = FilteredAlignment('HTR1A.ffn', 'HTR1A.faa', na_colCutoff=0.1)
    # qct = FilteredTree('HTR1A_P2N_na.aln')
    #
    # os.chdir('/work5/r2294/bin/Datasnakes/Orthologs/Align/QualityControl')
    # x = CodemlRun('HTR1A_P2N_na.aln', 'HTR1A_iqtree.nwk', 'paml.ctl')
    #
    # x.cml.run(verbose=True)