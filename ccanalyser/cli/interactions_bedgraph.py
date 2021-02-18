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
@click.option('--normalise', help='Normalised bedgraph (Correct for number of cis reads)', default=False, is_flag=True)
@click.option('--binsize', help='Binsize to use for converting bedgraph to evenly sized genomic bins', default=0)
@click.option('--gzip', help='Compress output using gzip', default=False, is_flag=True)
@click.option('--scale_factor', help='Scale factor to use for bedgraph normalisation', default=1e6, type=click.INT)
@click.option('--sparse/--dense', help='Scale factor to use for bedgraph normalisation', default=True)
def interactions_bedgraph(
    cooler_fn: str,
    capture_names: list = None,
    output_prefix: str = "",
    normalise: bool = False,
    binsize=0,
    gzip=True,
    scale_factor=1e6,
    sparse=True,
):

    """
    Converts cooler file to bedgraph.         

    Args:

    """

    if not capture_names:
        capture_names = [c.strip('/') for c in cooler.fileops.list_coolers(cooler_fn) if not 'resolutions' in c]


    for ii, capture_name in enumerate(capture_names):

        if binsize == 0:
            bedgraph = CoolerBedGraph(cooler_fn=f"{cooler_fn}::{capture_name}", sparse=sparse)

        elif ii == 0 and binsize > 0: # Only want to bin once and then re-use this for the rest
            binner = CoolerBinner(cooler_fn=f"{cooler_fn}::{capture_name}", binsize=binsize)
            bedgraph = CoolerBedGraphWindowed(cooler_fn=f"{cooler_fn}::{capture_name}", binner=binner, sparse=sparse)     
        
        elif binsize > 0:
            bedgraph = CoolerBedGraphWindowed(cooler_fn=f"{cooler_fn}::{capture_name}", binner=binner, sparse=sparse)

        
        bedgraph.to_file(
            f'{output_prefix}.{capture_name}.bedgraph{".gz" if gzip else ""}',
            normalise=normalise,
            scale_factor=scale_factor,
        )
