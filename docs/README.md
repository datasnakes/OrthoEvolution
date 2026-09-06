# OrthoEvolution documentation

The authored documentation lives in `docs/user_guide/`. Great Docs combines
these guides with the project README, package metadata, changelog, citation,
and an API reference generated from the Python docstrings.

## Prerequisites

- Python 3.11 through 3.14
- [Quarto](https://quarto.org/docs/get-started/)
- [`uv`](https://docs.astral.sh/uv/getting-started/installation/)
- A repository-local virtual environment

## Install the documentation tools

From the repository root:

```bash
uv venv --python 3.14 .venv
uv pip install --python .venv/bin/python -e ".[docs]"
```

## Build the site

```bash
.venv/bin/great-docs build --no-refresh
```

The generated site is written to `great-docs/_site/`. The entire
`great-docs/` directory is ephemeral and is not committed.

The `--no-refresh` option is intentional. OrthoEvolution contains bundled
BioSQL sources and Cookiecutter templates that make unconstrained package-wide
API discovery slow. The maintained API inventory in `great-docs.yml` keeps the
published reference explicit and reviewable.

## Preview the site

```bash
.venv/bin/great-docs preview
```

## Quality checks

```bash
.venv/bin/great-docs lint
.venv/bin/great-docs check-links --docs-only
.venv/bin/great-docs seo
```

GitHub Actions builds every documentation change. Pushes to `main` publish the
rendered site to GitHub Pages.
