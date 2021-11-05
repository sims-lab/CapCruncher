from typing import Iterable, List, Literal, Tuple
import dask
import pandas as pd

import click
import xopen
from capcruncher.cli.cli_alignments import cli
from capcruncher.utils import hash_column, load_json
import ujson
import os
import numpy as np
from dask.delayed import delayed
import dask.dataframe as dd
import shutil


def extract_start_and_end_slice_coords_pe(df):
    coords = df["coordinates"].str.extract(
        r"^chr(?P<chrom1>.*?):(?P<start>\d+).*\|chr(?P<chrom2>.*?):\d+-(?P<end>\d+)"
    )
    coords = coords["chrom1"].str.cat(coords.iloc[:, 1:])
    coords.name = "coordinates_pe"
    return coords


def identify_deduplicates_from_tsv(
    fragments: pd.DataFrame, read_type: Literal["flashed", "pe"], buffer: float = 1e6
):

    df_fragments = pd.read_csv(
        fragments,
        sep="\t",
        chunksize=buffer,
        usecols=["parent_read", "coordinates"],
    )

    unique_coordinates = set()
    duplicated_fragments = set()

    for ii, df in enumerate(fragments):

        # Shuffle to stop fragments at the end of the sample always being removed
        df = df.sample(frac=1)

        if read_type == "flashed":

            lst_id = hash_column(df["parent_read"])
            lst_coords = hash_column(df["coordinates"])

        else:
            # Extract chrom1 + start1 + chrom(last entry) + end(last entry)
            coords_df = df["coordinates"].str.extract(
                r"^chr(?P<chrom1>.*?):(?P<start>\d+).*\|chr(?P<chrom2>.*?):\d+-(?P<end>\d+)"
            )

            lst_id = hash_column(df["parent_read"])

            # chrom1+start1+chrom-1+end-1(hashed)
            lst_coords = hash_column(coords_df["chrom1"].str.cat(coords_df.iloc[:, 1:]))

        for id, coord in zip(lst_id, lst_coords):
            if not coord in unique_coordinates:
                unique_coordinates.add(coord)
            else:
                duplicated_fragments.add(id)

    return duplicated_fragments


def identify_duplicates_from_hdf5(
    fragments: list, read_type: Literal["flashed", "pe"]
):

    df_fragments_coords = dd.read_hdf(
        fragments, key=f"/fragments", columns=["id", "coordinates"]
    )

    if read_type == "flashed":

        return (
            df_fragments_coords.map_partitions(
                lambda df: df.assign(coordinates=hash_column(df["coordinates"]))
            )
            .shuffle(on="coordinates")
            .map_partitions(lambda df: df[df.duplicated(subset="coordinates")])
            ["id"]
        )

    elif read_type == "pe":

        return (
            df_fragments_coords.map_partitions(
                lambda df: df.join(extract_start_and_end_slice_coords_pe(df))
                [["id", "coordinates_pe"]]
                #.assign(coordinates_pe=lambda df: hash_column(df["coordinates_pe"]))
            )
            .shuffle(on="coordinates_pe")
            .map_partitions(lambda df: df[df.duplicated(subset="coordinates_pe")])
            ["id"]
        )


def identify_duplicates_from_parquet(
    fragments: list, viewpoint: str, read_type: Literal["flashed", "pe"]
):

    df_fragments_coords = dd.read_parquet(
        fragments,
        filters=[
            [
                ("viewpoint", "==", viewpoint),
            ]
        ],
        columns=["id", "coordinates"],
    )

    if read_type == "flashed":

        ids_duplicated = df_fragments_coords.shuffle(on="coordinates").map_partitions(
            lambda df: df[df.duplicated(subset="coordinates")]
        )["id"]

    elif read_type == "pe":

        ids_duplicated = (
            df_fragments_coords.map_partitions(
                lambda df: df.join(extract_start_and_end_slice_coords_pe(df))[
                    ["id", "coordinates_pe"]
                ]
            )
            .shuffle(on="coordinates_pe")
            .map_partitions(lambda df: df[df.duplicated(subset="coordinates_pe")])["id"]
        )

    return ids_duplicated


def remove_duplicates_from_tsv(
    slices: os.PathLike, output: os.PathLike, duplicated_ids: os.PathLike, buffer=1e6
):

    # Remove output if it exist as will need to append to file.
    if os.path.exists(output):
        os.unlink(output)

    df_slices = pd.read_csv(slices, sep="\t", chunksize=buffer)

    with xopen.xopen(duplicated_ids, "r") as r:
        ids_duplicated = {int(x) for x in ujson.load(r)}

    n_reads_total = 0
    n_reads_unique = 0

    # Iterate slices in chunks
    for ii, df in enumerate(df_slices):

        print(f"Processed {(ii + 1) * buffer} slices")

        n_reads_total += df["parent_read"].nunique()

        # Hash the parent_read column and remove any duplicated ids.
        df = (
            df.assign(parent_read_hashed=lambda df: hash_column(df["parent_read"]))
            .set_index("parent_read_hashed")
            .loc[lambda df: ~df.index.isin(ids_duplicated)]
        )

        n_reads_unique += df["parent_read"].nunique()

        # Append to file.
        df.reset_index(drop=True).to_csv(
            output, sep="\t", index=None, header=True if ii < 1 else False, mode="a"
        )

    return (n_reads_total, n_reads_unique)


def remove_duplicates_from_hdf5_single(
    slices: Iterable, duplicated_ids: os.PathLike, output: os.PathLike
):

    try:
        ser_duplicated_ids = pd.read_hdf(duplicated_ids, key="/duplicated_ids")
    except KeyError:
        ser_duplicated_ids = pd.Series(data=["NO_DATA"], name="duplicated_ids")

    with pd.HDFStore(slices, "r") as store:
        n_slices_total = store.get_storer(f"/slices").nrows
        df_slices_dedup = store.select(
            f"/slices", where="parent_id != ser_duplicated_ids"
        )
        n_slices_unique = df_slices_dedup.shape[0]

    df_slices_dedup.to_hdf(output, key=f"/slices", format="table")

    return (n_slices_total, n_slices_unique)

def remove_duplicates_from_hdf5_files (
    slices: Iterable, duplicated_ids: os.PathLike, output: os.PathLike
) -> Tuple[int, int]:

    try:
        ser_duplicated_ids = pd.read_hdf(duplicated_ids, key="/duplicated_ids")
    except KeyError:
        ser_duplicated_ids = pd.Series(data=["NO_DATA"], name="/duplicated_ids")

    n_slices_total = 0

    # Need to get total number of slices
    for slice_file in slices:
        with pd.HDFStore(slice_file, "r") as store:
            n_slices_total += store.get_storer("slices").nrows
    
    # Generate a dask dataframe
    dframes = [delayed(pd.read_hdf(fn, key="slices", where="parent_id != ser_duplicated_ids"))
               for fn in slices]

    ddf = dd.from_delayed(dframes, verify_meta=False)
    ddf.to_hdf(output, key="slices", format="table", data_columns=["viewpoint"], mode="w", min_itemsize={"slice_name": 75, "parent_read": 75})

    # Need to get final number of slices
    with pd.HDFStore(output, "r") as store:
        n_slices_unique = store.get_storer(f"slices").nrows

    return (n_slices_total, n_slices_unique)



def identify(
    fragments: os.PathLike,
    input_type: str = "parquet",
    output: os.PathLike = "duplicated_ids.json",
    viewpoint: str = "",
    buffer: int = 1e6,
    read_type: str = "flashed",
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

    if input_type == "tsv":
        if len(fragments) > 1:
            raise NotImplementedError("Currently just supports a single tsv input")
        else:
            duplicated_fragments = identify_deduplicates_from_tsv(
                fragments[0], read_type=read_type, buffer=buffer
            )
            with xopen.xopen(f"{output}", "w") as w:
                ujson.dump(dict.fromkeys(duplicated_fragments), w)

    elif input_type == "hdf5":

        outfile = f"{output.replace('.hdf5', '')}.hdf5"
        if os.path.exists(outfile):
            os.remove(outfile)

        duplicated_fragments = identify_duplicates_from_hdf5(
            fragments, read_type=read_type
        )

        duplicated_fragments.to_hdf(outfile, f"/duplicated_ids", min_itemsize={"id": 25})

    elif input_type == "parquet":

        import pyarrow.parquet as pq

        outfile = f"{output.replace('.parquet', '')}.parquet"

        if os.path.exists(outfile):
            shutil.rmtree(outfile)

        if viewpoint == "":
            ds = pq.ParquetDataset(fragments[0], use_legacy_dataset=False)
            viewpoints = pd.Series(ds.files).str.extract(r"viewpoint=(.*?)/").tolist()

        else:
            # Need a fake list to iterate
            viewpoints = [
                viewpoint,
            ]

        for viewpoint in viewpoints:
            duplicated_fragments = identify_duplicates_from_parquet(
                fragments, viewpoint=viewpoint, read_type=read_type
            )


def remove(
    slices: os.PathLike,
    duplicated_ids: os.PathLike,
    output: os.PathLike = "dedup.slices.tsv.gz",
    buffer: int = 5e6,
    sample_name: str = "",
    read_type: str = "",
    stats_prefix: os.PathLike = "",
    input_type: str = "hdf5",
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

    if input_type == "tsv":
        n_slices_total, n_slices_unique = remove_duplicates_from_tsv(
            slices, output, duplicated_ids, buffer=buffer
        )

    elif input_type == "hdf5":
        #slices = [slices,] if isinstance(slices, str) else slices 
        n_slices_total, n_slices_unique = remove_duplicates_from_hdf5_files(slices, duplicated_ids, output)

    # Prepare stats
    df_stats = pd.DataFrame()
    df_stats["stat_type"] = ["not-deduplicated", "deduplicated"]
    df_stats["stat"] = [n_slices_total, n_slices_unique]
    df_stats["sample"] = sample_name
    df_stats["read_type"] = read_type
    df_stats["read_number"] = 0
    df_stats["stage"] = "deduplicate_slices"
    df_stats.to_csv(f"{stats_prefix}.read.stats.csv", index=False)

    print(df_stats)
