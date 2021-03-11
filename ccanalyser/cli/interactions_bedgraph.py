#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

import cooler
import click
from ccanalyser.cli import cli
from ccanalyser.tools.pileup import CoolerBedGraph, CoolerBedGraphWindowed
from ccanalyser.tools.storage import CoolerBinner
import os


@cli.command()
@click.argument("cooler_fn")
@click.option(
    "-n",
    "--capture_names",
    help="Capture to extract and convert to bedgraph, if not provided will transform all.",
    multiple=True,
)
@click.option("-o", "--output_prefix", help="Output prefix for bedgraphs")
@click.option(
    "--normalise",
    help="Normalised bedgraph (Correct for number of cis reads)",
    default=False,
    is_flag=True,
)
@click.option(
    "--binsize",
    help="Binsize to use for converting bedgraph to evenly sized genomic bins",
    default=0,
)
@click.option("--gzip", help="Compress output using gzip", default=False, is_flag=True)
@click.option(
    "--scale_factor",
    help="Scale factor to use for bedgraph normalisation",
    default=1e6,
    type=click.INT,
)
@click.option(
    "--sparse/--dense",
    help="Produce bedgraph containing just positive bins (sparse) or all bins (dense)",
    default=True,
)
def interactions_bedgraph(
    cooler_fn: os.PathLike,
    capture_names: list = None,
    output_prefix: os.PathLike = "",
    normalise: bool = False,
    binsize: int = 0,
    gzip: bool = True,
    scale_factor: int = 1e6,
    sparse: bool = True,
):
    """Extracts a bedgraph for a specific capture probe/ all capture probes present in a cooler file.
       
       Bedgraph can be normalised for the number of cis interactions and/or binned into even genomic
       windows.

    Args:
     cooler_fn (os.PathLike): Path to hdf5 file containing cooler groups.
     capture_names (list, optional): Name of capture probe to extract.
                                     If None, will process all probes present in the file.
                                     Defaults to None.
     output_prefix (os.PathLike, optional): Output file prefix for bedgraph. Defaults to "".
     normalise (bool, optional): Normalise counts using the number of cis interactions. Defaults to False.
     binsize (int, optional): Genomic binsize to use for generating bedgraph. No binning performed if less than 0. Defaults to 0.
     gzip (bool, optional): Compress output bedgraph with gzip. Defaults to True.
     scale_factor (int, optional): Scaling factor for normalisation. Defaults to 1e6.
     sparse (bool, optional): Produce bedgraph containing just positive bins (True) or all bins (False). Defaults to True.
    """

    if (
        not capture_names
    ):  # If no probe names provided extract all restriction fragment cooler objects
        capture_names = [
            c.strip("/")
            for c in cooler.fileops.list_coolers(cooler_fn)
            if not "resolutions" in c
        ]

    for ii, capture_name in enumerate(capture_names):

        if binsize == 0:
            bedgraph = CoolerBedGraph(
                cooler_fn=f"{cooler_fn}::{capture_name}", sparse=sparse
            )

        elif (
            ii == 0 and binsize > 0
        ):  # Only want to bin once and then re-use this for the rest
            binner = CoolerBinner(
                cooler_fn=f"{cooler_fn}::{capture_name}", binsize=binsize
            )
            bedgraph = CoolerBedGraphWindowed(
                cooler_fn=f"{cooler_fn}::{capture_name}", binner=binner, sparse=sparse
            )

        elif binsize > 0:
            bedgraph = CoolerBedGraphWindowed(
                cooler_fn=f"{cooler_fn}::{capture_name}", binner=binner, sparse=sparse
            )

        bedgraph.to_file(
            f'{output_prefix}.{capture_name}.bedgraph{".gz" if gzip else ""}',
            normalise=normalise,
            scale_factor=scale_factor,
        )
