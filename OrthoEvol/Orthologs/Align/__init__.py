"""Alignment command line tool wrapper."""

from .orthoclustal import ClustalO
from .guidance2 import Guidance2Commandline
from .pal2nal import PAL2NALCommandline as Pal2NalCommandline
from .msa import MultipleSequenceAlignment

# Make this explicit, then they show up in the API docs
__all__ = ("ClustalO",
           "Guidance2Commandline",
           "Pal2NalCommandline",
           "MultipleSequenceAlignment",
           )
