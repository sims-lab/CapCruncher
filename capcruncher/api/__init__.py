"""CapCruncher public Python API.

Key entry points for programmatic use::

    from capcruncher.api import (
        digest_fastq,
        split_fastq,
        deduplicate_fastq,
        digest_genome,
        count_interactions,
        annotate_intervals,
    )

Imports are lazy: submodules (and their heavy deps) are only loaded on first
attribute access, not on ``import capcruncher.api``.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from capcruncher.api.fastq import deduplicate_fastq, digest_fastq, split_fastq
    from capcruncher.api.genome import digest_genome
    from capcruncher.api.interactions.count import count_interactions
    from capcruncher.api.intervals.annotate import annotate_intervals

__all__ = [
    "deduplicate_fastq",
    "digest_fastq",
    "split_fastq",
    "digest_genome",
    "count_interactions",
    "annotate_intervals",
]

_PUBLIC: dict[str, str] = {
    "deduplicate_fastq": "capcruncher.api.fastq",
    "digest_fastq": "capcruncher.api.fastq",
    "split_fastq": "capcruncher.api.fastq",
    "digest_genome": "capcruncher.api.genome",
    "count_interactions": "capcruncher.api.interactions.count",
    "annotate_intervals": "capcruncher.api.intervals.annotate",
}


def __getattr__(name: str) -> object:
    if name in _PUBLIC:
        module = importlib.import_module(_PUBLIC[name])
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
