# OrthoEvolution

[![CI](https://github.com/datasnakes/OrthoEvolution/actions/workflows/ci.yml/badge.svg)](https://github.com/datasnakes/OrthoEvolution/actions/workflows/ci.yml)
[![PyPI](https://badge.fury.io/py/OrthoEvol.svg)](https://pypi.org/project/OrthoEvol/)
[![Documentation](https://readthedocs.org/projects/orthoevolution/badge/?version=latest)](https://orthoevolution.readthedocs.io/en/latest/)
[![Coverage](https://codecov.io/gh/datasnakes/OrthoEvolution/branch/main/graph/badge.svg)](https://codecov.io/gh/datasnakes/OrthoEvolution)
[![DOI](https://zenodo.org/badge/88282824.svg)](https://doi.org/10.5281/zenodo.17796234)
[![Last Commit](https://badgen.net/github/last-commit/datasnakes/OrthoEvolution)](https://github.com/datasnakes/OrthoEvolution/commits/main)

OrthoEvolution is a Python package for reproducible comparative evolutionary
genetics, with a focus on ortholog inference, sequence analysis, and
phylogenetic workflows.

**Current version:** 1.0.0

## Table of Contents

- [Project Background](#project-background)
- [Core Capabilities](#core-capabilities)
- [Install & Setup](#install--setup)
- [Usage](#usage)
- [Documentation and Examples](#documentation-and-examples)
- [Testing](#testing)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Authors](#authors)

## Project Background

OrthoEvolution supports the inference and analysis of orthologous genes using
NCBI BLAST, multiple-sequence alignment strategies, and phylogenetic tools. It
organizes these steps into reusable workflows so researchers can manage large
comparative-genetics datasets and reproduce their analyses.

The package is organized around four major areas:

- `Orthologs` provides ortholog inference, alignment, and phylogenetic tools.
- `Manager` creates and coordinates repositories, projects, databases, and
  research datasets.
- `Tools` provides reusable utilities for data retrieval, parallel execution,
  logging, and cluster workloads.
- `Cookies` provides project and website templates.

For additional scientific context, see this
[related comparative-genetics paper](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2014.00283/full).

## Core Capabilities

- Infer candidate orthologs and generate post-BLAST reports.
- Retrieve NCBI datasets and preformatted BLAST databases.
- Prepare and filter nucleotide or protein sequence alignments.
- Support phylogenetic workflows involving PAML, PhyML, IQ-TREE, Phylip, and
  ETE.
- Create consistent directory structures for comparative-genetics projects.
- Configure local, parallel, PBS, and Slurm-oriented workloads.

Some workflows call external scientific programs or remote services. Install
the required BLAST, alignment, or phylogenetic software for the specific
workflow you intend to run.

## Install & Setup

OrthoEvolution supports Python 3.11 through 3.14. A virtual environment keeps its
dependencies separate from other Python projects.

Install [`uv`](https://docs.astral.sh/uv/getting-started/installation/) before
creating the environment.

### Install from PyPI

```bash
uv venv --python 3.14 .venv
uv pip install --python .venv/bin/python OrthoEvol
```

### Install from source

```bash
git clone https://github.com/datasnakes/OrthoEvolution.git
cd OrthoEvolution
uv venv --python 3.14 .venv
uv pip install --python .venv/bin/python .
```

### Install for development

```bash
git clone https://github.com/datasnakes/OrthoEvolution.git
cd OrthoEvolution
uv venv --python 3.14 .venv
uv pip install --python .venv/bin/python -e ".[test]"
```

## Usage

### Run a preconfigured local BLAST workflow

```python
from OrthoEvol.Orthologs.Blast import OrthoBlastN

gpcr_blastn = OrthoBlastN(
    project="orthology-gpcr",
    method=1,
    save_data=True,
    acc_file="gpcr.csv",
    copy_from_package=True,
)
gpcr_blastn.run()
```

This workflow requires a compatible local BLAST installation and database.

### Create a comparative-genetics project

```python
from OrthoEvol.Manager.management import ProjectManagement

project_manager = ProjectManagement(
    repo="test-repo",
    user=None,
    project="test-project",
    research=None,
    research_type="comparative_genetics",
    new_project=True,
)
```

### Download an NCBI BLAST database

```python
from pathlib import Path

from OrthoEvol.Tools.ftp import NcbiFTPClient

ncbi_ftp = NcbiFTPClient(email="researcher@example.org")
ncbi_ftp.getblastdb(
    database_name="refseq_rna",
    download_path=Path("databases"),
    v5=True,
)
```

NCBI database downloads require network access and can use substantial disk
space. Choose the destination and database deliberately before starting a
transfer.

## Documentation and Examples

- Read the
  [OrthoEvolution documentation](https://orthoevolution.readthedocs.io/en/latest/)
  for module and API details.
- Browse the [examples](examples/) for scripts, example data, and interface
  prototypes.
- Report problems or request enhancements through
  [GitHub Issues](https://github.com/datasnakes/OrthoEvolution/issues).

## Testing

Install the development dependencies and run the test suite through the active
virtual environment:

```bash
uv pip install --python .venv/bin/python -e ".[test]"
.venv/bin/python -m pytest tests/
```

The continuous-integration workflow runs the suite on Python 3.11, 3.12, 3.13,
and 3.14.

## Contributing

Contributions are welcome. Create a focused branch, include tests and
documentation where appropriate, and review the
[contributing guidelines](CONTRIBUTING.rst) before opening a pull request.

## Citation

If you use OrthoEvolution in research, please cite the software:

> Gilmore, R., & Hutchins, S. D. (2026). *OrthoEvolution* (Version 1.0.0)
> [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.17796234

OrthoEvolution builds on the work of the Biopython community. We thank its
developers and contributors and ask users to cite Biopython when it supports
their analyses:

> Cock, P. J. A., et al. (2009). Biopython: Freely available Python tools for
> computational molecular biology and bioinformatics. *Bioinformatics*,
> 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

## License

OrthoEvolution is distributed under the [MIT License](LICENSE).

## Authors

OrthoEvolution was created and is maintained by the Datasnakes:

- [Rob Gilmore](https://github.com/grabear)
- [Shaurita D. Hutchins](https://github.com/sdhutchins)
