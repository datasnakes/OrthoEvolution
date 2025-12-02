# docs

This directory is the source directory of our Sphinx documentation.

## Quick Start

See the main [docs/README.md](../README.md) for instructions on building documentation.

## Automation

Documentation is now automated:
- Version is automatically read from `setup.py`
- API docs are regenerated using `sphinx-apidoc`
- Use `python docs/update_docs.py` or `make -C docs all` to update and build
