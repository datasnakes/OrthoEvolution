from Datasnakes.Orthologs.Align.QualityControl.filter import FilteredAlignment
from Datasnakes.Orthologs.Phylogenetics.IQTree.best_tree import FilteredTree
from Datasnakes.Orthologs.Phylogenetics.PAML.codeml import CodemlRun
import os
import time
import subprocess
from shutil import copy
from pathlib import Path


class OrthologPipeline(object):

    def __init__(self, genes, qsub_template, worker_template, home=os.getcwd()):
        self.index = Path(home)  # index directory
        self.raw_data = self.index.parent / Path('raw_data')
        self.data = self.index.parent / Path('data')
        self.qsub_template = self.index / qsub_template
        self.worker_template = self.index / worker_template
        self.genes = genes

    def iterate(self):
        with open(self.genes, 'r') as gene_list:
            for gene in gene_list.readlines():
                command = self.batch_script_setup(str(self.qsub_template), str(self.worker_template), str(self.raw_data / Path(gene)), 'rgilmore',
                                        'rgilmore@gmail.com', '8GB', gene)
                self.submit(command)

    @staticmethod
    def batch_script_setup(self, qsub_file, python_file, raw_data_path, author, email, memory, gene):
        date_time = time.strftime("%H:%M:%S on %d/%m/%Y")
        qsub_command = "qsub -v PYTHONFILE=%s RAWDATA=%s AUTHOR=%s EMAIL=%s GBS=%s GENE=%s DATE=%s %s" % \
                       (python_file, raw_data_path, author, email, memory, gene, date_time, qsub_file)
        print('Command:  %s' % qsub_command)
        return qsub_command

    @staticmethod
    def submit(self, qsub_command):
        subprocess.check_call([qsub_command], stderr=subprocess.STDOUT, shell=True)
