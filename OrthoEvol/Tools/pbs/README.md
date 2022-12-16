# PBS Documentation

The Portable Batch System (or PBS for short) is a popular piece of computer software for Unix clusters.  It is used to
allocate hardware resources for computational tasks.  The PBS module in this package is used to connect with and partially
wrap some of the PBS command line tools for use in Python3 scripts.  Qstat and Qsub are currently the only tools that are
available in this module, but they provide plenty of functionality to get PBS jobs up and running through the Python3
interface.

Qsub can be used to create jobs on a cluster that contains the PBS software, while Qstat can be used to monitor or watch
existing jobs on the cluster.

## Qsub

Qsub is OrthoEvol's solution for submitting a PBS job using Python.  It primarily uses the _qsub <pbs_script>_ command to 
submit jobs.  This class also lets the end user utilize [string templating](https://docs.python.org/3.4/library/string.html#string.Template.)
for Python or PBS scripts.  Custom PBS scripts can also be created using the optional Qsub parameters.

#### BaseQsub class Example:

```python
from OrthoEvol.Tools.pbs import BaseQsub

# Initialize the job before submitting, including directory creation.
base_job = BaseQsub(job_name="GH_test", pbs_script="test.pbs")
# Copy the supplied pbs script to the new directory named with the <job_name> parameter
base_job.copy_supplied_script(supplied_script=base_job.supplied_pbs_script, new_script=base_job.pbs_script)
# Submit the PBS script
base_job.submit_pbs_script()
# Use the PBS job_id created by the PBS system
print(base_job.pbs_job_id)
```

#### Qsub class Example:

```python
from OrthoEvol.Tools.pbs import Qsub

# For Custom PBS scripts:
    # Initialize the job before submitting, including directory creation.
job = Qsub(python_script='test.py', job_name='GH_test',author='grabear', description='This is an example on GitHub.', 
           pbs_command_list=["pyenv activate pbs37"])
# Set up new directory with PBS/python scripts, and submit the job
job.submit_python_job()
# To resubmit the job use the rerun parameter
job.submit_python_job(rerun=True)

# For templated PBS and/or Python scripts:
    # If you have created files with string templating, then set up the attributes for python and PBS separately
...
py_attrs = {"pbs_wd": job.pbs_working_dir, "author":"grabear", "iterations": 4}
pbs_attrs = {"pbs_wd": job.pbs_working_dir, "author":"grabear", "walltime": "72:00:00"}
# Use them in the job submission call
job.submit_python_job(py_template_file="test.py", python_attributes=py_attrs, pbs_template_file="test.pbs", 
                      pbs_attributes=pbs_attrs)
# Use the PBS job_id created by the PBS system                    
print(job.pbs_job_id)
```

## Qstat

Qstat was developed in order to help monitor PBS jobs.  It's primary function is to collect data on a submitted job every
interval (which is specified by the user).  While we do have a proprietary means of parsing qstat output, more modern versions
of PBSpro have an option to parse qstat data into JSON format directly.

### Examples

#### BaseQstat

```python
from OrthoEvol.Tools.pbs import BaseQstat
from pprint import pprint
# Initialize qstat
job = BaseQstat(job_id="12345.sequoia")
# Run qstat
job.run_qstat(csv_flag=False)

# print data
pprint(job.qstat_dict)
```

#### Qstat

```python
from OrthoEvol.Tools.pbs import Qstat
from pprint import pprint
import pandas as pd

# Initialize qstat with an interval of 240 seconds
job = Qstat(job_id="12345.sequoia", wait_time=240)
# Watch the job and collect data for 10 intervals
job.watch(max_count=10)
# Print the plottable data
data = pd.DataFrame.from_csv(str(job.data_file))
data = data.to_dict()
pprint(data)

```
