# -*- coding: utf-8 -*-
import logging
import luigi
import os

from OrthoEvol.Tools.sge import SGEPipelineTask
from OrthoEvol.Orthologs.Blast import OrthoBlastN


# This is more pythonic with YAML loading
blast_config = {
  "taxon_file": None,
  "go_list": None,
  "post_blast": True,
  "template": None,
  "save_data": True,
  "copy_from_package": True,
  "MAF": 'test_blast.csv'
}


myblast = OrthoBlastN(proj_mana=None, project="sdh-test", project_path=os.getcwd(), **blast_config)


# TIP Works on linux
logger = logging.getLogger('luigi-interface')

SGEPipelineTask.shared_tmp_dir = os.getcwd()
SGEPipelineTask.parallel_env = None


class BlastPipelineTask(SGEPipelineTask):
    """Task for running a BLAST pipeline on SunGrid Engine.

    This task runs BLAST searches in parallel across multiple cluster nodes.
    It inherits from SGEPipelineTask to enable distributed execution.

    :param path: Working directory path for the pipeline.
    :type path: str
    :param accessions: List of accession numbers to process.
    :type accessions: str
    """

    path = luigi.Parameter()
    accessions = luigi.Parameter()

    def run(self):
        """Execute the BLAST pipeline task.

        Configures and runs BLAST for Homo sapiens using the configured
        blast_human settings. Override work() instead of run() for debugging.
        """
        myblast.blast_config(myblast.blast_human, 'Homo_sapiens', auto_start=True)

    def output(self):
        """Define the output target for this task.

        :return: Local target representing the project path output.
        :rtype: luigi.LocalTarget
        """
        return luigi.LocalTarget(path=os.path.join(self.path, myblast.project_path.as_posix()))


if __name__ == '__main__':
    path = os.getcwd()
    accessions = myblast.acc_list[1:]
    num_accs = len(accessions)
    tasks = [BlastPipelineTask(path=path,
                               accessions=str(accessions),
                               select=i+1) for i in range(num_accs)]
    luigi.build(tasks, local_scheduler=True, workers=num_accs)
