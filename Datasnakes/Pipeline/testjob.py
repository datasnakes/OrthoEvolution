import logging
import luigi
import os
from Datasnakes.Pipeline import SGEJobTask

# TIP Works on linux

logger = logging.getLogger('luigi-interface')

SGEJobTask.shared_tmp_dir = os.getcwd()
SGEJobTask.parallel_env = None

class MyJobTask(SGEJobTask):

    i = luigi.Parameter()

    def work(self):
        logger.info('Running test job...')
        with open(self.output().path, 'w') as f:
            f.write('This is a test job.')
        f.close()

    def output(self):
        return luigi.LocalTarget(path=os.path.join(os.getcwd(), 'testfile_' + str(self.i)))


if __name__ == '__main__':
    tasks = [MyJobTask(i=str(i), select=i+1) for i in range(3)]
    luigi.build(tasks, local_scheduler=True, workers=3)
