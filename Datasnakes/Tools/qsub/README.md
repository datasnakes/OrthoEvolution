qsub Documentation
-------------------------
Collection of tools for using PBS, a job scheduler for high-performance
computing environments. The command is usually `qsub <options>` on most systems.


Usage
-----
The main class under Genbank is `CreateJob`. Some functions are `import_temp`,
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


Tests
-----

Describe and show how to run the tests with code examples.

:exclamation: Notes
-------------------

Explain or list any notable information about the contents of this folder.

Other
-----
