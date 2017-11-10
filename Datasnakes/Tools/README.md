Tools Documentation
=====================
These tools were created by or modified our team to aid with the Orthologs package.

We've incorporated tools for bash with pybasher, qsub tools for use with pbs, a pandoc
script for converting docx files to markdown formats, multiprocessing in multiprocess, and
a ftp module that aids in downloading files from NCBI's ftp repository.


Examples
---------

### Download NCBI databases with our NCBI FTP Client
``` python
import Datasnakes.Tools.ftp import NcbiFTPClient
```

### Utilize multiprocessing to speed up your code
``` python
import Datasnakes.Tools
```

### Integrate logging in a simple and quick way
``` python
from Datasnakes.Tools.sge import SGEJob
```


#### More Information

Check the specific modules for more detailed readmes and examples of using the
tools with this package.