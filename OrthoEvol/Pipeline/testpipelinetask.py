import logging
import luigi
import os
from OrthoEvol.Tools.sge import SGEPipelineTask

# TIP Works on linux
logger = logging.getLogger('luigi-interface')

SGEPipelineTask.shared_tmp_dir = os.getcwd()
SGEPipelineTask.parallel_env = None


class TestPipelineTask(SGEPipelineTask):
    """Example pipeline task for testing SGE pipeline functionality.

    This is a simple example task that creates a test output file.
    It demonstrates the basic structure of a pipeline task.

    :param i: Task identifier number.
    :type i: str
    """

    i = luigi.Parameter()

    def work(self):
        """Execute the test pipeline task.

        Creates a test output file with a simple message.
        Use work() instead of run() for SGE tasks.
        """
        logger.info('Running test job...')
        with open(self.output().path, 'w') as f:
            f.write('This is a test job.')

    def output(self):
        """Define the output target for this task.

        :return: Local target representing the test output file.
        :rtype: luigi.LocalTarget
        """
        return luigi.LocalTarget(path=os.path.join(os.getcwd(), 'testjob_' + str(self.i) + '.txt'))


if __name__ == '__main__':
    tasks = [TestPipelineTask(i=str(i), select=i+1) for i in range(3)]
    luigi.build(tasks, local_scheduler=True, workers=3)
