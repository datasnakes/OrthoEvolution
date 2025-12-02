# Documentation

This directory contains the Sphinx documentation for OrthoEvolution.

## Quick Start

To update and build the documentation:

```bash
# From the project root
python docs/update_docs.py
```

This script will:
1. Automatically read the version from `setup.py`
2. Update Sphinx configuration
3. Convert README.md files to RST format (if pypandoc installed)
4. Regenerate API documentation from source code
5. Build the HTML documentation

### Options

```bash
# Skip README.md conversion (faster, if READMEs haven't changed)
python docs/update_docs.py --skip-readmes

# Only regenerate API docs, don't build
python docs/update_docs.py --skip-build
```

## Manual Steps

If you prefer to build manually:

### 1. Regenerate API Documentation

```bash
cd docs
rm -rf docs/source/modules/*.rst
sphinx-apidoc OrthoEvol/ -o docs/source/modules --separate --force
```

### 2. Build Documentation

```bash
cd docs/docs/source
make html
```

Or using sphinx-build directly:

```bash
sphinx-build -b html docs/source docs/_build
```

## Automation Features

### Version Management

The `conf.py` file automatically reads the version from `setup.py`, so you don't need to manually update version numbers when releasing.

### API Documentation

The `update_docs.py` script automatically:
- Removes old module documentation files
- Regenerates API docs using `sphinx-apidoc`
- Builds the final documentation

## Documentation Structure

- `docs/source/` - Sphinx source files
  - `conf.py` - Sphinx configuration (auto-reads version from setup.py)
  - `index.rst` - Main documentation index
  - `modules/` - Auto-generated API documentation
  - `cookies/`, `manager/`, `orthologs/`, `pipeline/`, `tools/` - Module-specific docs
  - `tutorial/` - Tutorial documentation

## Software Dependencies

- [Sphinx](http://www.sphinx-doc.org/) - Documentation generator
- [Pandoc](http://johnmacfarlane.net/pandoc/) - Optional, for converting markdown to RST

Install dependencies:

```bash
pip install sphinx
# Optional: pip install pypandoc
```

## ReadTheDocs Integration

The documentation is automatically built and hosted on [ReadTheDocs](http://orthoevolution.readthedocs.io/).

The `.readthedocs.yml` file in the project root configures the build process.

## Notes

- The `_static/` directory contains generated HTML files that are updated when docs are built
- README.md to RST conversion is integrated into `update_docs.py` (requires pypandoc)
- Module documentation is auto-generated, so manual edits to `modules/*.rst` files will be overwritten
