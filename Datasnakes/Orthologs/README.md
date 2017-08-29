Orthologs Documentation
-------------------------
This top level module includes submodules such as [Align](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Align/README.md) (for aligning multi fasta files),
[Phylogenetics](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Phylogenetics/README.md) (for analyzing multiple sequence alignments), [BioSQL]() (for database creation),
[Blast](https://github.com/datasnakes/Datasnakes-Scripts/tree/master/Datasnakes/Orthologs/Blast) (includes tools for using NCBI's blastn command line), and [Genbank](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Genbank/README.md).
(for tools to extract features from genbank files).

Usage
-----
These classes are optimized to be used together (very little work to do that),
but can also be used singularly.


#### Simple Example

This is a simple example of using some of the modules.

``` python
from Datasnakes.Orthologs import Phylogenetics as Phylo

```

Software Dependencies
----------------------
Ensure that the following software is installed and in your path:
1. Clustal omega
2. NCBI Standalone Blast
3. PAML
4. PhyML
5. Phylip

If you are a sudo user, you may use the script we've provided, [sudo-install.sh](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/sudo-install.sh).

Using `sudo-install.sh` on Debian/Ubuntu:

``` bash
# Change to the directory of the file.
cd
chmod +x sudo-install.sh
./sudo-install.sh

```

If you don't have sudo privileges, it'd be best to have your system admin
install the packages unless you are familiar with compiling packages from a source.