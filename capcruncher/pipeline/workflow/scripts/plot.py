# ruff: noqa: F821

import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger

import capcruncher.api.plotting as cp

logger.add(snakemake.log[0], format="{time} {level} {message}", level="INFO")

with logger.catch():
    can_group_tracks = (
        True if snakemake.params.design.groupby("condition").size().max() > 1 else False
    )

    # Set-up tracks
    tracks = []

    # Add scale bar
    tracks.append(cp.CCTrack(None, file_type="scale"))

    # Bigwig tracks
    if snakemake.input.bigwigs:
        if can_group_tracks:
            df_bw = pd.DataFrame(
                [pathlib.Path(p) for p in snakemake.input.bigwigs], columns=["fn"]
            )
            df_bw = df_bw.assign(
                samplename_and_vp=lambda df: df.fn.apply(lambda x: x.stem),
                samplename=lambda df: df.samplename_and_vp.str.replace(
                    f"_{snakemake.params.viewpoint}", ""
                ),
            ).merge(
                snakemake.params.design,
                left_on="samplename",
                right_on="sample",
                how="left",
            )

            for condition, df in df_bw.groupby("condition"):
                tracks.append(
                    cp.CCTrack(
                        df.fn.tolist(),
                        file_type="bigwig_summary",
                        title=condition,
                        min=0,
                        max="auto",
                    )
                )
                tracks.append(cp.CCTrack(None, file_type="spacer"))

        else:
            for bw in snakemake.input.bigwigs:
                bw_path = pathlib.Path(bw)
                tracks.append(
                    cp.CCTrack(
                        bw, file_type="bigwig", title=bw_path.stem, min=0, max="auto"
                    )
                )
                tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add subtractions if available
    if snakemake.input.subtractions:
        for sub in snakemake.input.subtractions:
            sub_path = pathlib.Path(sub)
            tracks.append(
                cp.CCTrack(sub, file_type="bigwig", title=sub_path.stem, min="auto")
            )
            tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add heatmaps if available
    if snakemake.input.heatmaps:

        if can_group_tracks:

            df_hm = pd.DataFrame(
                [pathlib.Path(p) for p in snakemake.input.heatmaps], columns=["fn"]
            )
            df_hm = df_hm.assign(
                samplename_and_vp=lambda df: df.fn.apply(lambda x: x.stem),
                samplename=lambda df: df.samplename_and_vp.str.replace(
                    f"_{snakemake.params.viewpoint}", ""
                ),
            ).merge(
                snakemake.params.design,
                left_on="samplename",
                right_on="sample",
                how="left",
            )

            for condition, df in df_hm.groupby("condition"):
                tracks.append(
                    cp.CCTrack(
                        df.fn.tolist(),
                        file_type="heatmap_summary",
                        title=condition,
                        binsize=snakemake.params.binsize,
                        viewpoint=snakemake.params.viewpoint,
                        style="triangular",
                    )
                )
                tracks.append(cp.CCTrack(None, file_type="spacer"))
        else:

            for hm in snakemake.input.heatmaps:
                hm_path = pathlib.Path(hm)
                tracks.append(
                    cp.CCTrack(
                        hm,
                        file_type="heatmap",
                        title=hm_path.stem,
                        binsize=snakemake.params.binsize,
                        viewpoint=snakemake.params.viewpoint,
                        style="triangular",
                    )
                )
                tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add genes if available
    if snakemake.params.genes:
        genes = snakemake.params.genes
        genes_path = pathlib.Path(genes)
        tracks.append(cp.CCTrack(genes, file_type="genes"))
        tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add X-axis
    tracks.append(cp.CCTrack(None, file_type="xaxis"))

    # Make figure and save
    fig = cp.CCFigure(tracks)
    fig.save(snakemake.params.coordinates, output=snakemake.output.fig)

    # Export template used to make figure
    fig.to_toml(snakemake.output.template)
