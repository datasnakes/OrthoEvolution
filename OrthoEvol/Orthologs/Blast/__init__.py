from .orthologs_blastn import OrthoBlastN
from .base_blastn import BaseBlastN
from .comparative_genetics import BaseComparativeGenetics, ComparativeGenetics
from OrthoEvol.Orthologs import OrthologsWarning
import warnings

# Make this explicit, then they show up in the API docs
__all__ = ("BaseComparativeGenetics",
           "ComparativeGenetics",
           "OrthoBlastN",
           "BaseBlastN"
           )

warnings.warn('Ensure that `BLASTDB` and `WINDOW_MASKER_PATH` are set as environment variables.',OrthologsWarning)