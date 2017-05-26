"""Genbank tools."""

from datasnakes.Orthologs.Genbank.archive.db2fasta import Database2Fasta as DB2FASTA

# Make this explicit, then they show up in the API docs
__all__ = ("DB2GBK",
           "DB2FASTA",
           )
