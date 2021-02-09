#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

import cooler
import click
from ccanalyser.cli import cli
from ccanalyser.tools.pileup import CoolerBedGraph, CoolerBedGraphWindowed
from ccanalyser.tools.storage import CoolerBinner

@cli.command()
@click.argument('cooler_fn')
@click.option('-n', '--capture_names', help='Capture to extract and convert to bedgraph, if not provided will transform all.', multiple=True)
@click.option('-o', '--output_prefix', help='Output prefix for bedgraphs')
@click.option('--normalise', help='Normalised bedgraph (Correct for number of cis reads)', default=False)
@click.option('--binsizes', help='Binsizes to use for converting bedgraph to evenly sized genomic bins', multiple=True)
@click.option('--gzip', help='Compress output using gzip', default=False)
@click.option('--scale_factor', help='Scale factor to use for bedgraph normalisation', default=1e6, type=click.INT)
def reporters_bedgraph(
    cooler_fn: str,
    capture_names: list = None,
    output_prefix: str = "",
    normalise: bool = False,
    binsize=0,
    gzip=True,
    normalise_scale=1e6,
):

    """
    Converts cooler file to bedgraph.         

    Args:

    """

    if not capture_names:
        capture_names = cooler.fileops.listcooler_fns(cooler_fn)

    if binsize > 0:
        # Set up binner object as don't want to run binning twice
        window = True
        binner = CoolerBinner(
            cooler_fn=f"{cooler_fn}::{capture_names[0]}",
            binsize=binsize,
            normalise=False,
        )

    for capture_name in capture_names:

        if window:
            bedgraph = CoolerBedGraphWindowed(
                cooler_fn=f"{cooler_fn}::{capture_name}", binner=binner
            )
        else:
            bedgraph = CoolerBedGraph(cooler_fn=f"{cooler_fn}::{capture_name}")

        bedgraph.to_file(
            f'{output_prefix}{capture_name}.bedgraph{".gz" if gzip else ""}',
            normalise=normalise,
            scale_factor=normalise_scale,
        )
