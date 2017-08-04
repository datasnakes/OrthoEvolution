Qsub Documentation
-------------------------
Collection of tools for using PBS, a job scheduler for high-performance
computing environments. The command is usually `qsub <options>` on most systems.


Usage
-----
The base class under qsub is `QsubUtils`. Some functions are `import_temp`,
which allows the user to import a preformatted template pbs script or python
script and use it in the pipeline if needed.

The class currently provides a template, `temp.pbs`, file to be modified and used
with the class using the `pbs_dict` function.

#### Code Examples

##### Submit 1 job

``` python

```

##### Submit multiple jobs

``` python

```

:exclamation: Notes
-------------------
Explain or list any notable information about the contents of this folder.
