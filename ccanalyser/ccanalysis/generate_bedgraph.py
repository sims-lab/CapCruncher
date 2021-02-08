#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

import cooler
from ccanalyser.ccanalysis.bedgraph import CoolerBedGraph, CoolerBedGraphWindowed
from ccanalyser.ccanalysis.storage import CoolerBinner


def main(
    cooler_file: str,
    capture_names: list = None,
    output_prefix: str = "",
    normalise: bool = False,
    binsize=0,
    gzip=True,
    normalise_scale=1e6,
):

    """
    Converts cooler file to bedgraph

    Args:

    """

    if not capture_names:
        capture_names = cooler.fileops.list_coolers(cooler_file)

    if binsize > 0:
        # Set up binner object as don't want to run binning twice
        window = True
        binner = CoolerBinner(
            cooler_fn=f"{cooler_file}::{capture_names[0]}",
            binsize=binsize,
            normalise=False,
        )

    for capture_name in capture_names:

        if window:
            bedgraph = CoolerBedGraphWindowed(
                cooler_fn=f"{cooler_file}::{capture_name}", binner=binner
            )
        else:
            bedgraph = CoolerBedGraph(cooler_fn=f"{cooler_file}::{capture_name}")

        bedgraph.to_file(
            f'{output_prefix}{capture_name}.bedgraph{".gz" if gzip else ""}',
            normalise=normalise,
            scale_factor=normalise_scale,
        )
