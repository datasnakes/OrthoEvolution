from .orthologs_blastn import OrthoBlastN
from .comparative_genetics_files import CompGenFiles
from .comparative_genetics_objects import CompGenObjects
from Datasnakes.Orthologs import OrthologsWarning
import warnings

# Make this explicit, then they show up in the API docs
__all__ = ("CompGenObjects",
           "CompGenFiles",
           "OrthoBlastN"
           )

warnings.warn('Ensure that `BLASTDB` and `WINDOW_MASKER_PATH` are set as environment variables.',OrthologsWarning)