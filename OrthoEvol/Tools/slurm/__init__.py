"""Submit Slurm jobs and inspect scheduler records."""

from .client import SlurmClient, SlurmCommandNotFoundError, SlurmJob

__all__ = ("SlurmClient", "SlurmCommandNotFoundError", "SlurmJob")
