# -*- coding: utf-8 -*-
import logging
import luigi
import os

from OrthoEvol.Tools.sge import SGEPipelineTask
from OrthoEvol.Orthologs.Blast import OrthoBlastN


# This is more pythonic with YAML loading
Blast_config = {
  "taxon_file": None,
  "go_list": None,
  "post_blast": True,
  "template": None,
  "save_data": True,
  "copy_from_package": True,
  "MAF": 'test_blast.csv'
}


myblast = OrthoBlastN(proj_mana=None, project="sdh-test", project_path=os.getcwd(), **Blast_config)


# TIP Works on linux
logger = logging.getLogger('luigi-interface')

SGEPipelineTask.shared_tmp_dir = os.getcwd()
SGEPipelineTask.parallel_env = None

class BlastPipelineTask(SGEPipelineTask):

    path = luigi.Parameter()
    accessions = luigi.Parameter()

    def run(self):  # Use work instead of run to DEBUG
        myblast.blast_config(myblast.blast_human, 'Homo_sapiens', auto_start=True)

    def output(self):
        return luigi.LocalTarget(path=os.path.join(self.path, myblast.project_path.as_posix()))


if __name__ == '__main__':
    path = os.getcwd()
    accessions = myblast.acc_list[1:]
    num_accs = len(accessions)
    tasks = [BlastPipelineTask(path=path,
                               accessions=str(accessions),
                               select=num_accs+1) for num_accs in range(num_accs)]
    luigi.build(tasks, local_scheduler=True, workers=num_accs)

