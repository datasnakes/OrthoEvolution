Multiprocess
=============
The multiprocess module and class, `MultiPro`, use python's native multiprocessing
module. Find more information [here](https://docs.python.org/3.6/library/multiprocessing.html).

In order to take advantage of using our supercomputer's processing power, we
looked into mpi and multiprocessing. Both were found to be useful.

This is an optional class in our pipeline, but if you're using AWS or Google's
supercomputing, then you may find it useful.

Usage
------

```python
from Datasnakes.Tools import Multiprocess

# Write a function that can be used for 1 item in a list of items
def blast2clustal(gene):

# Map your function to your list and run it using multiple processes
Multiprocess(n=5, blast2clustal(), geneslist)
```