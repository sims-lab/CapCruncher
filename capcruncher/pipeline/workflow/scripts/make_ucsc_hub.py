# ruff: noqa: F821

from __future__ import annotations

import pathlib
import re
import sys
from collections.abc import Iterable

from loguru import logger


SUMMARY_TRACK_PATTERN = re.compile(
    r"^(?P<sample>[^.]+)\.(?P<aggregation>[^.]+)-summary\.(?P<viewpoint>[^.]+)\.bigWig$"
)
COMPARISON_TRACK_PATTERN = re.compile(
    r"^(?P<sample>[^.]+-[^.]+)\.(?P<aggregation>[^.]+)-subtraction\.(?P<viewpoint>[^.]+)\.bigWig$"
)
REPLICATE_TRACK_PATTERN = re.compile(
    r"^(?P<sample>.+)_(?P<viewpoint>[^_]+)\.bigWig$"
)


def capcruncher_track_metadata(path: pathlib.Path) -> dict[str, str]:
    """Extract TrackNado metadata from CapCruncher result paths."""
    basename = path.name
    metadata: dict[str, str] = {}

    if path.suffix.lower() in {".bb", ".bigbed"}:
        return {
            "category": "Annotation",
            "normalisation": "viewpoints",
            "sample": "viewpoints",
            "viewpoint": "viewpoints",
            "aggregation": "viewpoints",
            "name": "viewpoint",
        }

    if path.suffix != ".bigWig":
        raise ValueError(f"Unsupported UCSC hub track type: {path}")

    summary_match = SUMMARY_TRACK_PATTERN.fullmatch(basename)
    if summary_match:
        metadata.update(summary_match.groupdict())
        metadata["category"] = "Aggregated"
        metadata["normalisation"] = "norm"
    else:
        comparison_match = COMPARISON_TRACK_PATTERN.fullmatch(basename)
        if comparison_match:
            metadata.update(comparison_match.groupdict())
            metadata["category"] = "Subtraction"
            metadata["normalisation"] = "norm"
        else:
            replicate_match = REPLICATE_TRACK_PATTERN.fullmatch(basename)
            if not replicate_match:
                raise ValueError(
                    "Could not parse CapCruncher track path. Expected one of: "
                    "<sample>_<viewpoint>.bigWig, "
                    "<sample>.<aggregation>-summary.<viewpoint>.bigWig, or "
                    "<sampleA>-<sampleB>.<aggregation>-subtraction.<viewpoint>.bigWig. "
                    f"Got: {path}"
                )
            metadata.update(replicate_match.groupdict())
            metadata["category"] = "Replicates"
            metadata["normalisation"] = path.parent.stem
            metadata["aggregation"] = "replicate"

    metadata["overlay"] = metadata["sample"]
    return metadata


def collect_track_paths(
    *,
    bigwigs: Iterable[str | pathlib.Path],
    bigwigs_summary: Iterable[str | pathlib.Path],
    bigwigs_comparison: Iterable[str | pathlib.Path],
    viewpoints: str | pathlib.Path,
) -> list[pathlib.Path]:
    """Return all files that should be added to the TrackNado hub builder."""
    return [
        *[pathlib.Path(path) for path in bigwigs],
        *[pathlib.Path(path) for path in bigwigs_summary],
        *[pathlib.Path(path) for path in bigwigs_comparison],
        pathlib.Path(viewpoints),
    ]


def build_track_metadata(
    *,
    bigwigs: Iterable[str | pathlib.Path],
    bigwigs_summary: Iterable[str | pathlib.Path],
    bigwigs_comparison: Iterable[str | pathlib.Path],
    viewpoints: str | pathlib.Path,
) -> list[dict[str, str]]:
    """Create the TrackNado metadata table for CapCruncher hub generation."""
    paths = collect_track_paths(
        bigwigs=bigwigs,
        bigwigs_summary=bigwigs_summary,
        bigwigs_comparison=bigwigs_comparison,
        viewpoints=viewpoints,
    )
    records = []
    for path in paths:
        metadata = capcruncher_track_metadata(path)
        records.append(
            {
                "fn": str(path),
                "basename": path.name,
                "ext": "bigBed"
                if path.suffix.lower() in {".bb", ".bigbed"}
                else "bigWig",
                **metadata,
            }
        )
    return records


def configure_logger(snakemake):
    """Route script logs to the Snakemake log file when one is available."""
    if not getattr(snakemake, "log", None):
        return

    logger.remove()
    logger.add(snakemake.log[0], format="{time} {level} {message}", level="INFO")
    logger.add(sys.stderr, format="{time} {level} {message}", level="ERROR")


def build_hub(
    *,
    tracks: Iterable[str | pathlib.Path],
    color_by: str,
    genome: str,
    hub_name: str,
    hub_email: str,
    outdir: str | pathlib.Path,
    report: str | pathlib.Path,
    custom_genome: bool | None = None,
    genome_twobit: str | pathlib.Path | None = None,
    genome_organism: str | None = None,
    genome_default_position: str | None = None,
):
    custom_genome_twobit: str | pathlib.Path = ""
    custom_genome_organism = genome
    if custom_genome:
        if not genome_twobit:
            raise ValueError(
                "Custom UCSC hub genomes require a genome twoBit file. "
                "Set genome.twobit in capcruncher_config.yml."
            )
        custom_genome_twobit = genome_twobit
        custom_genome_organism = genome_organism or genome

    import tracknado as tn

    builder = (
        tn.HubBuilder()
        .add_tracks([pathlib.Path(path) for path in tracks])
        .with_metadata_extractor(capcruncher_track_metadata)
        .group_by("category", "normalisation", as_supertrack=True)
        .group_by("sample", "viewpoint", "aggregation")
        .overlay_by("overlay")
        .color_by(color_by or "sample")
    )

    if custom_genome:
        builder = builder.with_custom_genome(
            name=genome,
            twobit_file=custom_genome_twobit,
            organism=custom_genome_organism,
            default_position=genome_default_position or "chr1:10000-20000",
        )

    return builder.build(
        name=hub_name,
        genome=genome,
        outdir=outdir,
        hub_email=hub_email,
        description_html=pathlib.Path(report),
    )


def main(snakemake):
    configure_logger(snakemake)
    logger.info("Getting data for UCSC hub tracks")
    tracks = collect_track_paths(
        bigwigs=snakemake.input.bigwigs,
        bigwigs_summary=snakemake.input.bigwigs_summary,
        bigwigs_comparison=snakemake.input.bigwigs_comparison,
        viewpoints=snakemake.input.viewpoints,
    )

    hub = build_hub(
        tracks=tracks,
        color_by=snakemake.params.color_by,
        genome=snakemake.params.genome,
        hub_name=snakemake.params.hub_name,
        hub_email=snakemake.params.hub_email,
        custom_genome=snakemake.params.custom_genome,
        genome_twobit=snakemake.params.genome_twobit,
        genome_organism=snakemake.params.genome_organism,
        genome_default_position=snakemake.params.genome_default_position,
        report=snakemake.input.report,
        outdir=snakemake.output[0],
    )
    hub.stage_hub()
    logger.info(f"Hub successfully generated in {snakemake.output[0]}")


if "snakemake" in globals():
    main(globals()["snakemake"])
