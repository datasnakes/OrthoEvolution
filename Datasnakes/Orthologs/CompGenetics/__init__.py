"""Comparative Genetics Analysis classes."""

from .comp_gen import CompGenObjects
from .ncbi_blast import CompGenFiles

# Make this explicit, then they show up in the API docs
__all__ = ("CompGenObjects",
           "CompGenFiles",
           )
