LogIt
=======
Use the LogIt class to make logging very simple. This short and sweet class
wraps around [logzero]() which allows color coded logging. We created our own
default logger with a default dateformat, logformat, and logging level (default
is debug). Setting up a logger is an easy 1 liner.

1. Import the LogIt class and create a variable. ex: `logit = LogIt()`
2. Create your logger. ex: `blastn = logit.default('blastn', 'blastn.log')`
2. Start logging. ex: `blastn.error('Your refseq accession was not found')`

Multiple loggers can exist for the same logfile and multiple loggers can be set
up for one script.

Example
---------

#### Use logging with ETE3PAML

```python
from Datasnakes.Tools import LogIt
from Datasnakes.Orthologs.Phylogenetics import ETE3PAML

# Set up your loggers
logit = LogIt()

# Log to one file
logfile = 'align2paml.log'

align, paml = logit.default('alignlog', logfile), logit.default('pamllog', logfile)

# Start logging
align.info('hi')
paml.info('muah')

# Shutdown the loggers and delete the logfile
logit.deletelog(logfile=logfile)
```