from Datasnakes.Manager.index.pipline import OrthologPipeline

x = OrthologPipeline(genes='genes.txt', qsub_template='pipeline.sh', worker_template='worker.py')

x.iterate()