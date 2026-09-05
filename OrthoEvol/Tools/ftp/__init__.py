"""FTP connection and NCBI bulk-download clients."""

from .baseftp import BaseFTPClient
from .ncbiftp import NcbiFTPClient


__all__ = ("BaseFTPClient", "NcbiFTPClient")
