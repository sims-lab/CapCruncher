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
    # Set-up tracks
    tracks = []

    # Add scale bar
    tracks.append(cp.CCTrack(None, file_type="scale"))

    # Bigwig tracks
    if snakemake.input.bigwigs:

        for bw in snakemake.input.bigwigs:
            bw_path = pathlib.Path(bw)
            tracks.append(cp.CCTrack(bw, file_type='bigwig', title=bw_path.stem))
            tracks.append(cp.CCTrack(None, file_type='spacer'))

    # Add subtractions if available
    if snakemake.input.subtractions:
        for sub in snakemake.input.subtractions:
            sub_path = pathlib.Path(sub)
            tracks.append(
                cp.CCTrack(sub, file_type='bigwig', title=sub_path.stem, min="auto")
            )
            tracks.append(cp.CCTrack(None, file_type='spacer'))

    # Add heatmaps if available
    if snakemake.input.heatmaps:
        for hm in snakemake.input.heatmaps:
            hm_path = pathlib.Path(hm)
            tracks.append(
                cp.CCTrack(
                    hm,
                    file_type='heatmap',
                    title=hm_path.stem,
                    binsize=snakemake.params.binsize,
                    viewpoint=snakemake.params.viewpoint,
                    style="triangular",
                )
            )
            tracks.append(cp.CCTrack(None, file_type='spacer'))

    # Add genes if available
    if snakemake.params.genes:
        genes = snakemake.params.genes
        genes_path = pathlib.Path(genes)
        tracks.append(cp.CCTrack(genes, file_type='genes'))
        tracks.append(cp.CCTrack(None, file_type='spacer'))

    # Add X-axis
    tracks.append(cp.CCTrack(None, file_type='xaxis'))

    # Make figure and save
    fig = cp.CCFigure(tracks)
    fig.save(snakemake.params.coordinates, output=snakemake.output.fig)

    # Export template used to make figure
    fig.to_toml(snakemake.output.template)
