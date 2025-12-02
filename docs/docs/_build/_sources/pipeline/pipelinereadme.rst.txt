Pipeline Documentation
======================

The Pipeline module is designed to provide the user with easily callable
and command line usable pipelines that allow orthology inference to be
completed in a parallel fashion.

This module uses Luigi and SunGrid Engine (SGE) to distribute
computational tasks across cluster nodes. Tasks are designed to run on
clusters that use **pbspro or SunGrid Engine**.

Overview
--------

The Pipeline module provides pre-configured pipeline tasks that can be
executed in parallel on cluster computing systems. Currently, the module
includes:

- **BlastPipelineTask**: Runs BLAST searches in parallel across multiple
  nodes
- **TestPipelineTask**: Example task for testing pipeline functionality

Examples
--------

Running a Test Pipeline
~~~~~~~~~~~~~~~~~~~~~~~

The ``TestPipelineTask`` is a simple example that demonstrates how to
create and run a pipeline task:

.. code:: python

   import logging
   import luigi
   import os
   from OrthoEvol.Tools.sge import SGEPipelineTask
   from OrthoEvol.Pipeline.testpipelinetask import TestPipelineTask

   # Configure SGE settings
   SGEPipelineTask.shared_tmp_dir = os.getcwd()
   SGEPipelineTask.parallel_env = None

   # Create and run test tasks
   tasks = [TestPipelineTask(i=str(i), select=i+1) for i in range(3)]
   luigi.build(tasks, local_scheduler=True, workers=3)

Running a BLAST Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~

The ``BlastPipelineTask`` runs BLAST searches in parallel:

.. code:: python

   import logging
   import luigi
   import os
   from OrthoEvol.Tools.sge import SGEPipelineTask
   from OrthoEvol.Pipeline.blastpipeline import BlastPipelineTask
   from OrthoEvol.Orthologs.Blast import OrthoBlastN

   # Configure BLAST settings
   blast_config = {
       "taxon_file": None,
       "go_list": None,
       "post_blast": True,
       "template": None,
       "save_data": True,
       "copy_from_package": True,
       "MAF": 'test_blast.csv'
   }

   # Initialize BLAST instance
   myblast = OrthoBlastN(
       proj_mana=None,
       project="sdh-test",
       project_path=os.getcwd(),
       **blast_config
   )

   # Configure SGE settings
   logger = logging.getLogger('luigi-interface')
   SGEPipelineTask.shared_tmp_dir = os.getcwd()
   SGEPipelineTask.parallel_env = None

   # Create and run BLAST tasks
   path = os.getcwd()
   accessions = myblast.acc_list[1:]
   num_accs = len(accessions)
   tasks = [
       BlastPipelineTask(
           path=path,
           accessions=str(accessions),
           select=i+1
       ) for i in range(num_accs)
   ]
   luigi.build(tasks, local_scheduler=True, workers=num_accs)

Task Parameters
---------------

All pipeline tasks inherit from ``SGEPipelineTask`` and support the
following parameters:

- **select**: Number of CPUs (slots) to allocate for the task (default:
  3)
- **shared_tmp_dir**: Shared drive accessible from all cluster nodes
  (default: ‘/home’)
- **parallel_env**: SGE parallel environment name (default: ‘orte’)
- **job_name**: Explicit job name for qsub
- **run_locally**: Run locally instead of on the cluster (default:
  False)

Software Dependencies
---------------------

- **Luigi**: Workflow management library
- **SunGrid Engine (SGE)**: Job scheduler for cluster computing
- **pbspro**: Alternative job scheduler (version 14.1.0 or higher)

Notes
-----

- Tasks should override the ``work()`` method instead of ``run()`` for
  SGE execution
- Use ``local_scheduler=True`` for local testing and debugging
- Set ``workers`` parameter to the number of parallel tasks you want to
  run
- Ensure Luigi is installed on all cluster nodes
