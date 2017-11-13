Orthologs Documentation
========================
This top level module includes submodules such as [Align](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Align/README.md) (for aligning multi fasta files),
[Phylogenetics](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Phylogenetics/README.md) (for analyzing multiple sequence alignments), [BioSQL]() (for database creation),
[Blast](https://github.com/datasnakes/Datasnakes-Scripts/tree/master/Datasnakes/Orthologs/Blast) (includes tools for using NCBI's blastn command line), and [Genbank](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/Genbank/README.md).
(for tools to extract features from genbank files).

## Usage & Examples
These classes are optimized to be used together (very little work to do that),
but can also be used as standalone classes/methods.

This is a simple example of using all of the `Orthologs` submodules together.


``` python
from Datasnakes.Orthologs.Blast import OrthoBlastN
from Datasnakes.Orthologs.Align import ClustalO
from Datasnakes.Orthologs.Phlogenetics import ETE3PAML

```

## ❗ Software Dependencies
Ensure that the following software is installed and in your path:
- Clustal omega
- NCBI Blast+ 2.6.0 or greater
- PAML
- PhyML
- Phylip
- IQTREE
- Mafft
- Prank
- Clustalw
- Guidance2
- Pal2Nal

If you are a sudo user, you may use the script we've provided, [install.sh](https://github.com/datasnakes/Datasnakes-Scripts/blob/master/Datasnakes/Orthologs/install.sh).

## Using `install.sh` on Debian/Ubuntu:

``` bash
# Change to the directory of the file.
cd
chmod +x install.sh
./sudo-install.sh

```
