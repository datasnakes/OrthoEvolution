# PBS Documentation

The Portable Batch System (or PBS for short) is a popular piece of computer software for Unix clusters.  It is used to
allocate hardware resources for computational tasks.  The PBS module in this packaged is used to connect with and partially
wrap some of the PBS command line tools for use in Python3 scripts.  Qstat and Qsub are currently the only tools that are
available in this module, but they provide plenty of functionality to get PBS jobs up and running through the Python3
interface.

Qsub can be used to create jobs on a cluster that contains the PBS software, while Qstat can be used to monitor or watch
existing jobs on the cluster.

## Qsub

Qsub is OrthoEvol's solution for submitting a PBS job using Python.  It primarily uses the _qsub <pbs_script>_ command to 
submit jobs.  This class also lets the end user utilize [string templating](https://docs.python.org/3.4/library/string.html#string.Template.)
for Python or PBS scripts.  Custom PBS scripts can also be created using the optional Qsub parameters.



## Qstat