"""NCBI BLASTn tool"""

from .blastn import BLASTn

# Make this explicit, then they show up in the API docs
__all__ = ("BLASTn",
)