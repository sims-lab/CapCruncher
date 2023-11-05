"""
Aim: Check that there is only one restriction fragment per viewpoint.
"""

import pandas as pd
import numpy as np
import pyranges as pr
import polars as pl
import pathlib
from loguru import logger
from capcruncher.api.annotate import BedIntersector


with logger.catch():
    logger.info("Checking that there is only one restriction fragment per viewpoint")

    bins = snakemake.input.bins
    viewpoints = snakemake.input.viewpoints
    
    gr = BedIntersector(viewpoints, bins, "restriction_fragments", 0.51).get_intersection(method="count")
    multiple_fragments = (gr.df["restriction_fragments"] > 1).sum()
    has_multiple_frags = multiple_fragments > 0


    if (
        has_multiple_frags
        and not snakemake.params.ignore_multiple_bins_per_viewpoint
    ):

        raise ValueError(
            f"""The following viewpoints overlap multiple restriction fragments:\n{gr}\n"""
        )

    else:
        pathlib.Path(snakemake.output.sentinel).touch()
