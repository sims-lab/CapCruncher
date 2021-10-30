#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from typing import Literal, Union

import cooler
from capcruncher.tools.pileup import CoolerBedGraph, CoolerBedGraphWindowed
from capcruncher.tools.storage import CoolerBinner


def bedgraph(
    cooler_fn: os.PathLike,
    viewpoint_names: list = None,
    output_prefix: os.PathLike = "",
    normalisation: Literal["raw", "n_cis", "region"] = "raw",
    normalisation_regions: os.PathLike = None,
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
     viewpoint_names (list, optional): Name of viewpoints to extract. 
                                       If None, will process all probes present in the file.
                                       Defaults to None.
     output_prefix (os.PathLike, optional): Output file prefix for bedgraph. Defaults to "".
     normalisation (bool, optional): Normalise counts using the number of cis interactions. Defaults to False.
     binsize (int, optional): Genomic binsize to use for generating bedgraph. No binning performed if less than 0. Defaults to 0.
     gzip (bool, optional): Compress output bedgraph with gzip. Defaults to True.
     scale_factor (int, optional): Scaling factor for normalisation. Defaults to 1e6.
     sparse (bool, optional): Produce bedgraph containing just positive bins (True) or all bins (False). Defaults to True.
    """

    viewpoint_names = viewpoint_names or [
        v.strip("/")
        for v in cooler.fileops.list_coolers(cooler_fn)
        if not "resolutions" in v
    ]

    bin_bedgraph = True if binsize > 0 else False

    for ii, viewpoint_name in enumerate(viewpoint_names):

        if not bin_bedgraph:
            bedgraph = CoolerBedGraph(
                cooler_fn=f"{cooler_fn}::{viewpoint_name}", sparse=sparse
            )

        elif ii == 0 and bin_bedgraph:
            # Only want to bin once and then re-use this for the rest
            binner = CoolerBinner(
                cooler_group=f"{cooler_fn}::{viewpoint_name}", binsize=binsize
            )
            bedgraph = CoolerBedGraphWindowed(
                cooler_fn=f"{cooler_fn}::{viewpoint_name}", binner=binner, sparse=sparse
            )

        else:
            bedgraph = CoolerBedGraphWindowed(
                cooler_fn=f"{cooler_fn}::{viewpoint_name}", binner=binner, sparse=sparse
            )

        if normalisation in ["n_cis", "region"]:
            bedgraph.normalise_bedgraph(
                scale_factor=scale_factor,
                method=normalisation,
                region=normalisation_regions,
            )

        bedgraph.to_file(
            f'{output_prefix}.{viewpoint_name}.bedgraph{".gz" if gzip else ""}',
        )
