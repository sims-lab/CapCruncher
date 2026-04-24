# ruff: noqa: F821

from __future__ import annotations

import pathlib
from collections.abc import Iterable

import pandas as pd
from loguru import logger


def _paths_to_frame(paths: Iterable[str | pathlib.Path], category: str) -> pd.DataFrame:
    return pd.DataFrame({"fn": [pathlib.Path(path) for path in paths]}).assign(
        basename=lambda df: df["fn"].map(lambda path: path.name),
        category=category,
    )


def _replicate_tracks(bigwigs: Iterable[str | pathlib.Path]) -> pd.DataFrame:
    df = _paths_to_frame(bigwigs, "Replicates")
    if df.empty:
        return df

    df["normalisation"] = df["fn"].map(lambda path: path.parent.stem)
    df[["sample", "viewpoint"]] = df["basename"].str.extract(
        r"(?P<sample>.*)_(?P<viewpoint>.*?).bigWig"
    )
    df["aggregation"] = "replicate"
    return df


def _summary_tracks(bigwigs: Iterable[str | pathlib.Path]) -> pd.DataFrame:
    df = _paths_to_frame(bigwigs, "Aggregated")
    if df.empty:
        return df

    df["normalisation"] = "norm"
    df[["sample", "aggregation", "viewpoint"]] = df["basename"].str.extract(
        r"(?P<sample>.*)\.(?P<aggregation>.*)-summary\.(?P<viewpoint>.*).bigWig"
    )
    return df


def _comparison_tracks(bigwigs: Iterable[str | pathlib.Path]) -> pd.DataFrame:
    df = _paths_to_frame(bigwigs, "Subtraction")
    if df.empty:
        return df

    df["normalisation"] = "norm"
    df[["sample", "aggregation", "viewpoint"]] = df["basename"].str.extract(
        r"(?P<sample>.*?)\.(?P<aggregation>.*?)-subtraction\.(?P<viewpoint>.*?).bigWig"
    )
    return df


def build_track_metadata(
    *,
    bigwigs: Iterable[str | pathlib.Path],
    bigwigs_summary: Iterable[str | pathlib.Path],
    bigwigs_comparison: Iterable[str | pathlib.Path],
    viewpoints: str | pathlib.Path,
) -> pd.DataFrame:
    """Create the TrackNado metadata table for CapCruncher hub generation."""
    tracks = pd.concat(
        [
            _replicate_tracks(bigwigs),
            _summary_tracks(bigwigs_summary),
            _comparison_tracks(bigwigs_comparison),
        ],
        ignore_index=True,
    )

    if not tracks.empty:
        tracks["overlay"] = tracks["sample"]
        tracks["ext"] = "bigWig"

    viewpoint_track = pd.DataFrame(
        [
            {
                "fn": pathlib.Path(viewpoints),
                "basename": pathlib.Path(viewpoints).name,
                "category": "Annotation",
                "normalisation": "viewpoints",
                "sample": "viewpoints",
                "viewpoint": "viewpoints",
                "aggregation": "viewpoints",
                "overlay": pd.NA,
                "ext": "bigBed",
                "name": "viewpoint",
            }
        ]
    )

    df = pd.concat([tracks, viewpoint_track], ignore_index=True)
    df["fn"] = df["fn"].map(str)
    return df


def build_hub(
    *,
    track_metadata: pd.DataFrame,
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
    import tracknado as tn

    builder = tn.HubBuilder().add_tracks_from_df(track_metadata)
    for track in getattr(builder, "tracks", []):
        if track.track_type != "bigWig":
            track.metadata.pop("overlay", None)

    builder = (
        builder.group_by("category", "normalisation", as_supertrack=True)
        .group_by("sample", "viewpoint", "aggregation")
        .overlay_by("overlay")
        .color_by(color_by)
    )

    if custom_genome:
        builder = builder.with_custom_genome(
            name=genome,
            twobit_file=genome_twobit,
            organism=genome_organism,
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
    logger.info("Getting data for UCSC hub tracks")
    track_metadata = build_track_metadata(
        bigwigs=snakemake.input.bigwigs,
        bigwigs_summary=snakemake.input.bigwigs_summary,
        bigwigs_comparison=snakemake.input.bigwigs_comparison,
        viewpoints=snakemake.input.viewpoints,
    )

    hub = build_hub(
        track_metadata=track_metadata,
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


if "snakemake" in globals():
    main(snakemake)
