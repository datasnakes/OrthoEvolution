# Slurm Tools

The Slurm tools submit an existing batch script and retrieve one scheduler
snapshot. They do not generate site-specific resource requests or poll the
scheduler.

```python
from pathlib import Path

from OrthoEvol.Tools.slurm import SlurmClient


client = SlurmClient()
job_id = client.submit(Path("ortholog-analysis.sh"))
active_jobs = client.active_jobs()
completed_job = client.job_history(job_id)
```

`submit()` uses `sbatch --parsable`. `active_jobs()` requests an explicit,
pipe-delimited format from `squeue`. `job_history()` uses `sacct --parsable2`
so it does not depend on optional JSON support.

Each method makes one command call. Callers control when another scheduler
snapshot is requested.
