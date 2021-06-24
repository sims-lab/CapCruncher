#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cooler
from capcruncher.tools.pileup import CoolerBedGraph, CoolerBedGraphWindowed
from capcruncher.tools.storage import CoolerBinner
import os

def bedgraph(
    cooler_fn: os.PathLike,
    capture_names: list = None,
    output_prefix: os.PathLike = "",
    normalise: bool = False,
    binsize: int = 0,
    gzip: bool = True,
    scale_factor: int = 1e6,
    sparse: bool = True,
):
    """
    Extracts reporters from a capture experiment and generates a bedgraph file.

    Identifies reporters for a single probe (if a probe name is supplied) or all capture 
    probes present in a capture experiment HDF5 file. 
    
    The bedgraph generated can be normalised by the number of cis interactions for
    inter experiment comparisons and/or binned into even genomic windows.

    \f
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
