Orthologs Documentation
-------------------------
This top level module includes submodules such as Align (for aligning multi fasta files),
Phylogenetics (for analyzing multiple sequence alignments), BioSQL (for database creation),
Blast (includes tools for using NCBI's blastn command line), and Genbank
(for tools to extract features from genbank files).

Usage
-----
These classes are optimized to be used together (very little work to do that),
but can also be used singularly.


#### Simple Example

This is a simple example of using some of the modules.

``` python
from Orthologs import Phylogenetics as Phylo

```

Software Dependencies
----------------------
Ensure that the following software is installed and in your `PATH`:
    1. Clustal omega
    2. NCBI Standalone Blast
    3. PAML
    4. PhyML
    5. Phylip

If you are a sudo user, you may use the script we've provided, `sudo-install.sh`.

Using `sudo-install.sh` on Debian/Ubuntu:
```bash
# Change to the directory of the file.
1. chmod +x sudo-install.sh
2. ./sudo-install.sh
```

If you don't have sudo privileges, it'd be best to have your system admin
install the packages unless you are familiar with compiling packages from a source.