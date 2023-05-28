# ruff: noqa: F821

import os
import pandas as pd
import itertools
import numpy as np
import pathlib
import re
from loguru import logger
import trackhub
import tracknado

# Single bigwigs
df_bw = pd.DataFrame(
    [pathlib.Path(p) for p in snakemake.input.bigwigs],
    columns=["fn"],
)

df_bw["basename"] = df_bw["fn"].apply(lambda p: p.name)
df_bw["normalisation"] = df_bw["fn"].apply(lambda p: p.parent)
df_bw[["sample", "viewpoint"]] = df_bw["basename"].str.extract(
    "(?P<sample>.*?)_(?P<viewpoint>.*?).bigWig"
)
df_bw["category"] = "Replicates"

# Summarised bigwigs
df_bw_summary = pd.DataFrame(
    [pathlib.Path(p) for p in snakemake.input.bigwigs_summary],
    columns=["fn"],
)
df_bw_summary["basename"] = df_bw_summary["fn"].apply(lambda p: p.name)
df_bw_summary["normalisation"] = "norm"
df_bw_summary[["sample", "aggregation", "viewpoint"]] = df_bw_summary[
    "basename"
].str.extract("(?P<sample>.*?)\.(?P<method>.*?)(?P<viewpoint>.*?).bigWig")
df_bw_summary["category"] = "Aggregated"

# Compared bigwigs
df_bw_compared = pd.DataFrame(
    [pathlib.Path(p) for p in snakemake.input.bigwigs_comparison],
    columns=["fn"],
)
df_bw_compared["basename"] = df_bw_compared["fn"].apply(lambda p: p.name)
df_bw_compared["normalisation"] = "norm"
df_bw_compared[["sample", "aggregation", "viewpoint"]] = df_bw_compared[
    "basename"
].str.extract("(.*?)\.(.*?)-subtraction\.(.*?).bigWig")
df_bw_compared["category"] = "Subtraction"

# Combine dataframes
df = pd.concat([df_bw, df_bw_summary, df_bw_compared], axis=0)

# Create hub design
design = tracknado.TrackDesign.from_design(
    df,
    color_by=snakemake.params.color_by,
    subgroup_by=["sample", "viewpoint", "aggregation"],
    supergroup_by=[
        "category",
    ],
    overlay_by=[
        "sample",
    ],
)

hub = tracknado.HubGenerator(
    track_design=design,
    genome=snakemake.params.genome,
    hub_name=snakemake.params.hub_name,
    description_html=pathlib.Path(snakemake.input.report),
    hub_email=snakemake.params.hub_email,
    custom_genome=snakemake.params.custom_genome,
    genome_twobit=snakemake.params.genome_twobit,
    genome_organism=snakemake.params.genome_organism,
    genome_default_position=snakemake.params.genome_default_position,
    outdir=snakemake.output[0],
)

hub.trackdb.add_tracks(
    trackhub.Track(
        name="viewpoint",
        tracktype="bigBed",
        source=snakemake.input.viewpoints,
        visibility="dense",
        color="0,0,0",
        autoScale="off",
        maxHeightPixels="100:50:8",
        shortLabel="Viewpoint",
        longLabel="Viewpoint",
    )
)

hub.stage_hub()
