#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import os
import subprocess
import tempfile
from typing import Literal, Union

import cooler
from capcruncher.tools.pileup import CoolerBedGraph, CoolerBedGraphWindowed
from capcruncher.tools.storage import CoolerBinner


def pileup(
    uri: os.PathLike,
    viewpoint_names: list = None,
    output_prefix: os.PathLike = "",
    format: Literal["bedgraph", "bigwig"] = "bedgraph",
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
     uri (os.PathLike): Path to hdf5 file containing cooler groups.
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

    logging.info(f"Performing pileup for {viewpoint_names}")
    viewpoint_names = viewpoint_names or [
        v.strip("/") for v in cooler.fileops.list_coolers(uri) if not "resolutions" in v
    ]

    bin_bedgraph = True if binsize > 0 else False

    for ii, viewpoint_name in enumerate(viewpoint_names):

        cooler_group = f"{uri}::{viewpoint_name}"

        if bin_bedgraph:
            cooler_group = f"{cooler_group}/resolutions/{binsize}"

        try:
            cooler.fileops.is_cooler(cooler_group)
        except Exception as e:
            logging.info(f"Exception {e} occured while looking for: {viewpoint_name}")
            raise (f"Cannot find {viewpoint_name} in cooler file")

        bedgraph = CoolerBedGraph(uri=cooler_group, sparse=sparse).extract_bedgraph(
            normalisation=normalisation,
            region=normalisation_regions,
            scale_factor=scale_factor,
        )

        logging.info(f"Generated bedgraph for {viewpoint_name}")

        if format == "bedgraph":

            bedgraph.to_csv(
                f'{output_prefix}.{viewpoint_name}.bedgraph{".gz" if gzip else ""}',
                sep="\t",
                header=False,
                index=False,
            )
        
        elif format == "bigwig":
            
            clr = cooler.Cooler(cooler_group)

            with tempfile.NamedTemporaryFile() as chromsizes_tmp:
                with tempfile.NamedTemporaryFile() as bedgraph_tmp:
                    clr.chromsizes.to_csv(chromsizes_tmp, sep="\t", header=False)
                    bedgraph.to_csv(bedgraph_tmp, sep="\t", index=False, header=False)

                    result = subprocess.run(["bedGraphToBigWig", bedgraph_tmp.name, chromsizes_tmp.name, f"{output_prefix}.{viewpoint_name}.bigWig"])




            

