import logging
import luigi
import os
from OrthoEvol.Tools.sge import SGEPipelineTask

# TIP Works on linux
logger = logging.getLogger('luigi-interface')

SGEPipelineTask.shared_tmp_dir = os.getcwd()
SGEPipelineTask.parallel_env = None

class TestPipelineTask(SGEPipelineTask):

    i = luigi.Parameter()

    def work(self):  # Use work instead of run to DEBUG
        logger.info('Running test job...')
        with open(self.output().path, 'w') as f:
            f.write('This is a test job.')
            f.close()

    def output(self):
        return luigi.LocalTarget(path=os.path.join(os.getcwd(), 'testjob_' + str(self.i) + '.txt'))


if __name__ == '__main__':
    tasks = [TestPipelineTask(i=str(i), select=i+1) for i in range(3)]
    luigi.build(tasks, local_scheduler=True, workers=3)
