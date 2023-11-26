# ruff: noqa: F821

import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from loguru import logger

import capcruncher.api.plotting as cp

logger.add(open(snakemake.log[0], "w"))

with logger.catch():
    logger.info("Checking if we can group tracks by condition")
    can_group_tracks = (
        True if snakemake.params.design.groupby("condition").size().max() > 1 else False
    )

    logger.info("Setting up tracks")
    # Set-up tracks
    tracks = []

    # Add scale bar
    tracks.append(cp.CCTrack(None, file_type="scale"))

    # Bigwig tracks
    if snakemake.input.bigwigs:
        logger.info("Adding bigwig tracks")
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
                logger.info(f"Added {condition} bigwig track")
                tracks.append(cp.CCTrack(None, file_type="spacer"))

        else:
            for bw in snakemake.input.bigwigs:
                bw_path = pathlib.Path(bw)
                tracks.append(
                    cp.CCTrack(
                        bw, file_type="bigwig", title=bw_path.stem, min=0, max="auto"
                    )
                )
                logger.info(f"Added {bw_path.stem} bigwig track")
                tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add subtractions if available
    if snakemake.input.subtractions:
        logger.info("Adding subtraction tracks")
        for sub in snakemake.input.subtractions:
            sub_path = pathlib.Path(sub)
            logger.info(f"Adding {sub_path.stem} subtraction track")
            tracks.append(
                cp.CCTrack(sub, file_type="bigwig", title=sub_path.stem, min="auto")
            )
            tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add heatmaps if available
    if snakemake.input.heatmaps:
        logger.info("Adding heatmaps")
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
                logger.info(f"Adding {condition} heatmap track")
                tracks.append(
                    cp.CCTrack(
                        df.fn.tolist(),
                        file_type="heatmap_summary",
                        title=condition,
                        binsize=snakemake.params.binsize,
                        viewpoint=snakemake.params.viewpoint,
                        style="triangular",
                        normalization=snakemake.params.normalization_method,
                    )
                )
                tracks.append(cp.CCTrack(None, file_type="spacer"))
        else:
            for hm in snakemake.input.heatmaps:
                hm_path = pathlib.Path(hm)
                logger.info(f"Adding {hm_path.stem} heatmap track")
                tracks.append(
                    cp.CCTrack(
                        hm,
                        file_type="heatmap",
                        title=hm_path.stem,
                        binsize=snakemake.params.binsize,
                        viewpoint=snakemake.params.viewpoint,
                        style="triangular",
                        normalization=snakemake.params.normalization_method,
                    )
                )
                tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add genes if available
    if snakemake.params.genes:
        logger.info("Adding genes track")
        genes = snakemake.params.genes
        genes_path = pathlib.Path(genes)
        tracks.append(cp.CCTrack(genes, file_type="genes"))
        tracks.append(cp.CCTrack(None, file_type="spacer"))

    # Add X-axis
    tracks.append(cp.CCTrack(None, file_type="xaxis"))

    # Make figure and save
    logger.info("Making figure")
    fig = cp.CCFigure(tracks)

    logger.info(f"Saving figure to: {snakemake.output.fig}")
    fig.save(snakemake.params.coordinates, output=snakemake.output.fig)

    # Export template used to make figure
    logger.info(f"Exporting template to {snakemake.output.template}")
    fig.to_toml(snakemake.output.template)
