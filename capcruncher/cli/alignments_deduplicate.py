import logging
from typing import Iterable, List, Literal, Tuple
import dask
import pandas as pd
import xopen
import ujson
import os
import numpy as np
import dask.distributed

from capcruncher.tools.deduplicate import (
    identify_coordinate_duplicates_from_tsv,
    identify_coordinate_duplicates_from_hdf5,
    identify_coordinate_duplicates_from_parquet,
    remove_duplicates_from_parquet,
    remove_duplicates_from_tsv,
    remove_duplicates_from_hdf5,
    read_duplicated_ids,
)
from capcruncher.utils import get_file_type


def identify(
    fragments: tuple,
    file_type: str = "auto",
    output: os.PathLike = "duplicated_ids.json",
    viewpoint: str = "",
    buffer: int = 1e6,
    read_type: str = "flashed",
    n_cores: int = 1,
    memory_limit: str = "1G",
):
    """
    Identifies aligned fragments with duplicate coordinates.

    Parses a tsv file containing filtered aligned fragments and generates a dictionary containing
    the hashed parental read id and hashed genomic coordinates of all slices. Duplicated fragments
    are implicitly removed if they share the same genomic coordinate hash.

    For non-combined reads (pe) a genomic coordinate hash is generated from the start of the first slice and
    the end of the last slice. This is due to the decreased confidence in the quality of the centre of the fragment.
    The coordinate hash for combined reads (flashed) is generated directly from the fragment coordinates. Only
    fragments with the exact coordinates and slice order will be considered to be duplicates.

    Identified duplicate fragments are output in json format to be used by the "remove" subcommand.

    \f
    Args:
     fragments_fn (os.PathLike): Input fragments.tsv file to process.
     output (os.PathLike, optional): Output path to output duplicated parental read ids. Defaults to "duplicated_ids.json".
     buffer (int, optional): Number of fragments to process in memory. Defaults to 1e6.
     read_type (str, optional): Process combined(flashed) or non-combined reads (pe).
                                Due to the low confidence in the quaility of pe reads, duplicates are identified by
                                removing any fragments with matching start and end coordinates.
                                Defaults to "flashed".
    """

    input_file_type = get_file_type(fragments[0]) if file_type == "auto" else file_type
    output_file_type = get_file_type(output)

    if input_file_type in ["hdf5", "parquet"]:
        # Will use dask for these
        cluster = dask.distributed.LocalCluster(
            n_workers=n_cores,
            threads_per_worker=1,
            dashboard_address=None,
            processes=True,
            scheduler_port=0,
            local_directory=os.environ.get("TMPDIR", "/tmp/"),
            memory_limit=memory_limit,
        )

    with dask.distributed.Client(cluster) as client:

        # Extract duplicated fragment ids
        if input_file_type == "tsv":
            if len(fragments) > 1:
                raise NotImplementedError("Currently just supports a single tsv input")
            else:
                duplicated_fragments = identify_coordinate_duplicates_from_tsv(
                    fragments[0], read_type=read_type, buffer=buffer
                )
        elif input_file_type == "hdf5":
            outfile = f"{output.replace('.hdf5', '')}.hdf5"
            if os.path.exists(outfile):
                os.remove(outfile)

            duplicated_fragments = identify_coordinate_duplicates_from_hdf5(
                fragments, read_type=read_type
            )

        elif input_file_type == "parquet":
            duplicated_fragments = identify_coordinate_duplicates_from_parquet(
                fragments, read_type=read_type
            )

        else:
            raise ValueError(f"Input file type {file_type} not supported")

        # Output
        if output_file_type == "json":
            with xopen.xopen(output, "w") as w:
                ujson.dump(dict.fromkeys(duplicated_fragments), w)

        elif output_file_type == "hdf5":

            duplicated_fragments.to_hdf(
                outfile, f"/duplicated_ids", min_itemsize={"id": 25}, index=False
            )

        elif output_file_type == "pickle":
            duplicated_fragments = duplicated_fragments.compute()
            duplicated_fragments.to_pickle(output)


def remove(
    slices: os.PathLike,
    duplicated_ids: os.PathLike,
    output: os.PathLike = "dedup.slices.tsv.gz",
    buffer: int = 5e6,
    sample_name: str = "",
    read_type: str = "",
    stats_prefix: os.PathLike = "",
    file_type: str = "hdf5",
    n_cores: int = 1,
    memory_limit: str = "1G",
):
    """
    Removes duplicated aligned fragments.

    Parses a tsv file containing aligned read slices and outputs only slices from unique fragments.
    Duplicated parental read id determined by the "identify" subcommand are located within the
    slices tsv file and removed.

    Outputs statistics for the number of unique slices and the number of duplicate slices identified.

    \f
    Args:
     slices_fn (os.PathLike): Input slices.tsv file.
     duplicated_ids (os.PathLike): Duplicated parental read ids in json format.
     output (os.PathLike, optional): Output file path for deduplicated slices. Defaults to "dedup.slices.tsv.gz".
     buffer (int, optional): Number of slices to process in memory. Defaults to 1e6.
     sample_name (str, optional): Name of sample being processed e.g. DOX-treated_1 used for statistics. Defaults to "".
     read_type (str, optional): Process combined(flashed) or non-combined reads (pe) used for statistics. Defaults to "".
     stats_prefix (os.PathLike, optional): Output path for deduplication statistics. Defaults to "".
    """

    input_file_type = (
        get_file_type(slices if isinstance(slices, str) else slices[0])
        if file_type == "auto"
        else file_type
    )
    output_file_type = get_file_type(output)

    duplicates = read_duplicated_ids(duplicated_ids)

    if input_file_type in ["hdf5", "parquet"]:
        # Will use dask for these
        cluster = dask.distributed.LocalCluster(
            n_workers=n_cores,
            threads_per_worker=1,
            dashboard_address=None,
            processes=True,
            scheduler_port=0,
            local_directory=os.environ.get("TMPDIR", "/tmp/"),
            memory_limit=memory_limit,
        )
    with dask.distributed.Client(cluster) as client:

        if input_file_type == "tsv":

            if not output_file_type == "tsv":
                raise NotImplementedError("Currently only tsv -> tsv output supported")

            n_reads_total, n_reads_unique = remove_duplicates_from_tsv(
                slices, output, duplicates, buffer=buffer
            )

        elif input_file_type == "hdf5":

            if not output_file_type == "hdf5":
                raise NotImplementedError("Currently only hdf5 -> hdf5 output supported")

            n_reads_total, n_reads_unique = remove_duplicates_from_hdf5(
                slices, duplicates, output
            )

        elif input_file_type == "parquet":
            if not output_file_type == "parquet":
                raise NotImplementedError(
                    "Currently only parquet -> parquet output supported"
                )

            n_reads_total, n_reads_unique = remove_duplicates_from_parquet(
                slices, duplicates, output
            )

    # Prepare stats
    df_stats = pd.DataFrame()
    df_stats["stat_type"] = ["not-deduplicated", "deduplicated"]
    df_stats["stat"] = [n_reads_total, n_reads_unique]
    df_stats["sample"] = sample_name
    df_stats["read_type"] = read_type
    df_stats["read_number"] = 0
    df_stats["stage"] = "deduplicate_slices"
    df_stats.to_csv(f"{stats_prefix}.read.stats.csv", index=False)

    return df_stats
