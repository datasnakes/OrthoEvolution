parallel
=============
The parellel module is home to the `Multiprocess` class which uses python's native multiprocessing
module. Find more information [here](https://docs.python.org/3.6/library/multiprocessing.html). It will
soon be home to [MPI (Message Passing Interface)](http://mpi4py.readthedocs.io/en/stable/) which is also a form of parallel computing.

In order to take advantage of using our supercomputer's processing power, we
looked into mpi and multiprocessing. Both were found to be useful.

This is an optional class in our pipeline, but if you're using AWS or Google's
supercomputing, then you may find it useful unless you're interested in or
using clustering or SGE (Sun Grid Engine). We have a sge module for that.

Examples
------

#### Random Example

```python
from Datasnakes.Tools import Multiprocess


def printwords(word):
    print(word)


words = ['bae', 'luh', 'cuh']

if __name__ == '__main__':
    mp = Multiprocess()
    mp.map2function(printwords, words)
```
