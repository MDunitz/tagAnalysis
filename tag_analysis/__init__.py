"""tag_analysis: importable 16S/18S amplicon (DADA2) processing pipeline.

Configure a run with RunConfig and call process_16s / process_18s.
"""

from .config import RunConfig, PrimerSet, PRIMERS_16S, PRIMERS_18S
from .pipelines import process, process_16s, process_18s
from . import constants

__version__ = "0.1.0"

__all__ = [
    "RunConfig",
    "PrimerSet",
    "PRIMERS_16S",
    "PRIMERS_18S",
    "process",
    "process_16s",
    "process_18s",
    "constants",
    "__version__",
]
