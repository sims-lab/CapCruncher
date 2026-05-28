"""CapCruncher public Python API.

Subpackages expose their own namespaced symbols::

    import capcruncher.api.alignments as alignments
    import capcruncher.api.filtering as filtering
    import capcruncher.api.interactions as interactions
    import capcruncher.api.intervals as intervals

Common entry points are also available directly::

    from capcruncher.api import (
        # fastq
        digest_fastq,
        split_fastq,
        deduplicate_fastq,
        FastqDigestOptions,
        FastqSplitOptions,
        FastqDeduplicationOptions,
        # genome
        digest_genome,
        # alignments
        annotate_alignments,
        filter_alignments,
        bam_to_parquet,
        parse_bam,
        AlignmentAnnotateOptions,
        AlignmentFilterOptions,
        # filtering
        CCSliceFilter,
        TriCSliceFilter,
        TiledCSliceFilter,
        FilterPipeline,
        FilterStepRegistry,
        # interactions
        count_interactions,
        deduplicate_interactions,
        cooler_to_bedgraph,
        pileup,
        differential,
        get_differential_interactions,
        summarise_reporter_viewpoints,
        write_countable_reporters,
        # intervals
        annotate_intervals,
        increase_cis_slice_priority,
        remove_duplicates_from_bed,
    )

Imports are lazy: submodules (and their heavy deps) are only loaded on first
attribute access, not on ``import capcruncher.api``.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # fastq
    # alignments
    from capcruncher.api.alignments.annotate import AlignmentAnnotateOptions
    from capcruncher.api.alignments.annotate import annotate as annotate_alignments
    from capcruncher.api.alignments.filter import (
        AlignmentFilterOptions,
        merge_annotations,
    )
    from capcruncher.api.alignments.filter import filter as filter_alignments
    from capcruncher.api.alignments.io import bam_to_parquet, parse_bam
    from capcruncher.api.fastq import (
        FastqDeduplicationOptions,
        FastqDigestOptions,
        FastqSplitOptions,
        deduplicate_fastq,
        digest_fastq,
        split_fastq,
    )

    # filtering
    from capcruncher.api.filtering.pipeline import (
        CCSliceFilter,
        FilterPipeline,
        SliceFilter,
        TiledCSliceFilter,
        TriCSliceFilter,
    )
    from capcruncher.api.filtering.steps import FilterStepName, FilterStepRegistry

    # genome
    from capcruncher.api.genome import digest_genome

    # interactions
    from capcruncher.api.interactions.bedgraph import (
        CCBedgraph,
        CoolerBedGraph,
        CoolerBedGraphWindowed,
        cooler_to_bedgraph,
    )
    from capcruncher.api.interactions.compare import concat as concat_interactions
    from capcruncher.api.interactions.compare import summarise as summarise_interactions
    from capcruncher.api.interactions.count import (
        InteractionCountOptions,
        count_interactions,
    )
    from capcruncher.api.interactions.deduplicate import (
        deduplicate as deduplicate_interactions,
    )
    from capcruncher.api.interactions.differential import (
        differential,
        get_differential_interactions,
    )
    from capcruncher.api.interactions.pileup import PileupOptions, pileup
    from capcruncher.api.interactions.reporters import (
        ReporterViewpointSummary,
        summarise_reporter_viewpoints,
        valid_viewpoint_names,
        write_countable_reporters,
    )

    # intervals
    from capcruncher.api.intervals.annotate import (
        annotate_intervals,
        increase_cis_slice_priority,
        remove_duplicates_from_bed,
    )

    # statistics models
    from capcruncher.api.statistics import (
        AlignmentDeduplicationStats,
        CisOrTransStats,
        DigestionStats,
        FastqDeduplicationStatistics,
        FastqTrimmingStatistics,
        SliceFilterStats,
        SliceFilterStatsList,
    )

__all__ = [
    # fastq
    "FastqDeduplicationOptions",
    "FastqDigestOptions",
    "FastqSplitOptions",
    "deduplicate_fastq",
    "digest_fastq",
    "split_fastq",
    # genome
    "digest_genome",
    # alignments
    "AlignmentAnnotateOptions",
    "AlignmentFilterOptions",
    "annotate_alignments",
    "bam_to_parquet",
    "filter_alignments",
    "merge_annotations",
    "parse_bam",
    # filtering
    "CCSliceFilter",
    "FilterPipeline",
    "FilterStepName",
    "FilterStepRegistry",
    "SliceFilter",
    "TiledCSliceFilter",
    "TriCSliceFilter",
    # interactions
    "CCBedgraph",
    "CoolerBedGraph",
    "CoolerBedGraphWindowed",
    "InteractionCountOptions",
    "PileupOptions",
    "ReporterViewpointSummary",
    "concat_interactions",
    "cooler_to_bedgraph",
    "count_interactions",
    "deduplicate_interactions",
    "differential",
    "get_differential_interactions",
    "pileup",
    "summarise_interactions",
    "summarise_reporter_viewpoints",
    "valid_viewpoint_names",
    "write_countable_reporters",
    # intervals
    "annotate_intervals",
    "increase_cis_slice_priority",
    "remove_duplicates_from_bed",
    # statistics
    "AlignmentDeduplicationStats",
    "CisOrTransStats",
    "DigestionStats",
    "FastqDeduplicationStatistics",
    "FastqTrimmingStatistics",
    "SliceFilterStats",
    "SliceFilterStatsList",
]

_PUBLIC: dict[str, str] = {
    # fastq
    "FastqDeduplicationOptions": "capcruncher.api.fastq",
    "FastqDigestOptions": "capcruncher.api.fastq",
    "FastqSplitOptions": "capcruncher.api.fastq",
    "deduplicate_fastq": "capcruncher.api.fastq",
    "digest_fastq": "capcruncher.api.fastq",
    "split_fastq": "capcruncher.api.fastq",
    # genome
    "digest_genome": "capcruncher.api.genome",
    # alignments
    "AlignmentAnnotateOptions": "capcruncher.api.alignments.annotate",
    "AlignmentFilterOptions": "capcruncher.api.alignments.filter",
    "bam_to_parquet": "capcruncher.api.alignments.io",
    "merge_annotations": "capcruncher.api.alignments.filter",
    "parse_bam": "capcruncher.api.alignments.io",
    # filtering
    "CCSliceFilter": "capcruncher.api.filtering.pipeline",
    "FilterPipeline": "capcruncher.api.filtering.pipeline",
    "FilterStepName": "capcruncher.api.filtering.steps",
    "FilterStepRegistry": "capcruncher.api.filtering.steps",
    "SliceFilter": "capcruncher.api.filtering.pipeline",
    "TiledCSliceFilter": "capcruncher.api.filtering.pipeline",
    "TriCSliceFilter": "capcruncher.api.filtering.pipeline",
    # interactions
    "CCBedgraph": "capcruncher.api.interactions.bedgraph",
    "CoolerBedGraph": "capcruncher.api.interactions.bedgraph",
    "CoolerBedGraphWindowed": "capcruncher.api.interactions.bedgraph",
    "InteractionCountOptions": "capcruncher.api.interactions.count",
    "PileupOptions": "capcruncher.api.interactions.pileup",
    "ReporterViewpointSummary": "capcruncher.api.interactions.reporters",
    "cooler_to_bedgraph": "capcruncher.api.interactions.bedgraph",
    "count_interactions": "capcruncher.api.interactions.count",
    "differential": "capcruncher.api.interactions.differential",
    "get_differential_interactions": "capcruncher.api.interactions.differential",
    "pileup": "capcruncher.api.interactions.pileup",
    "summarise_reporter_viewpoints": "capcruncher.api.interactions.reporters",
    "valid_viewpoint_names": "capcruncher.api.interactions.reporters",
    "write_countable_reporters": "capcruncher.api.interactions.reporters",
    # intervals
    "annotate_intervals": "capcruncher.api.intervals.annotate",
    "increase_cis_slice_priority": "capcruncher.api.intervals.annotate",
    "remove_duplicates_from_bed": "capcruncher.api.intervals.annotate",
    # statistics
    "AlignmentDeduplicationStats": "capcruncher.api.statistics",
    "CisOrTransStats": "capcruncher.api.statistics",
    "DigestionStats": "capcruncher.api.statistics",
    "FastqDeduplicationStatistics": "capcruncher.api.statistics",
    "FastqTrimmingStatistics": "capcruncher.api.statistics",
    "SliceFilterStats": "capcruncher.api.statistics",
    "SliceFilterStatsList": "capcruncher.api.statistics",
}

# Aliased names that differ from their source symbol
_ALIASES: dict[str, tuple[str, str]] = {
    "annotate_alignments": ("capcruncher.api.alignments.annotate", "annotate"),
    "filter_alignments": ("capcruncher.api.alignments.filter", "filter"),
    "concat_interactions": ("capcruncher.api.interactions.compare", "concat"),
    "deduplicate_interactions": (
        "capcruncher.api.interactions.deduplicate",
        "deduplicate",
    ),
    "summarise_interactions": ("capcruncher.api.interactions.compare", "summarise"),
}


def __getattr__(name: str) -> object:
    if name in _PUBLIC:
        module = importlib.import_module(_PUBLIC[name])
        return getattr(module, name)
    if name in _ALIASES:
        mod_path, attr = _ALIASES[name]
        module = importlib.import_module(mod_path)
        return getattr(module, attr)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
