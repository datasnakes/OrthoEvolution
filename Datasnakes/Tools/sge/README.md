sge Documentation
-------------------------
Collection of tools for using PBS, a job scheduler for high-performance
computing environments on SGE. The command is usually `qsub <options>` on most systems.

Usage
-----
The base class under qsub is `QsubUtils`. Some functions are `import_temp`,
which allows the user to import a preformatted template pbs script or python
script and use it in the pipeline if needed.

The class currently provides a template, `temp.pbs`, file to be modified and used
when submitting a job.

#### Code Examples

##### Submit 1 job

``` python

```

##### Submit multiple jobs

``` python

```

##### Get Job Info

``` python

```


:exclamation: Notes
-------------------
<br>

Thanks
-------------------
Thanks to [@jfeala](https://github.com/jfeala) for his work on Luigi's SGEJobTask.

