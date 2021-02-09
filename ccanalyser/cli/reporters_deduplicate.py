import pandas as pd

# import dask
# from dask.distributed import Client, LocalCluster
# import dask.array as da
# import dask.dataframe as dd
# import random

import click
import xopen
from ccanalyser.cli import cli
from ccanalyser.utils import hash_column, load_json
import numpy as np
import ujson
import os


@cli.group()
def reporters_deduplicate():
    """Identifies duplicate aligned fragments and removes them"""


@reporters_deduplicate.command()
@click.argument("fragments_fn")
@click.option(
    "-o",
    "--output",
    help="File to output fragments with duplicated coordinates",
    default="duplicated_ids.json.gz",
)
@click.option(
    "--buffer",
    help="Number of fragments to process at one time",
    default=1e6,
    type=click.INT,
)
def identify(fragments_fn, output="duplicated_ids.json", buffer=1e6):

    # TODO: Make the variable names cleaner and break this up a bit

    df_fragments = pd.read_csv(
        fragments_fn,
        sep="\t",
        chunksize=buffer,
        usecols=["parent_read", "pe", "coordinates"],
    )

    coordinates_flashed = dict()
    coordinates_pe = dict()
    parent_read_all = set()

    for df in df_fragments:
        flashed = df.query('pe == "flashed"')
        pe = df.query('pe == "pe"')

        # Deal with flashed coords
        # -> {coord(hashed): id(hashed)}
        flashed_dict = dict(
            zip(
                hash_column(flashed["coordinates"]), hash_column(flashed["parent_read"])
            )
        )

        # Extract chrom1 + start1 + chrom-1 + end-1
        pe_coords = pe["coordinates"].str.extract(
            r"^chr(?P<chrom1>[\d|X|Y|M]+):(?P<start>\d+).*\|chr(?P<chrom2>[\d|X|Y|M]+):\d+-(?P<end>\d+)"
        )

        # -> {chrom1+start1+chrom-1+end-1(hashed): id(hashed)}
        pe_dict = dict(
            zip(
                hash_column(pe_coords["chrom1"].str.cat(pe_coords.iloc[:, 1:])),
                hash_column(pe["parent_read"]),
            )
        )

        coordinates_flashed.update(flashed_dict)
        coordinates_pe.update(pe_dict)

        id_set = {x for x in [*flashed_dict.values(), *pe_dict.values()]}
        parent_read_all.update(id_set)

    # Identify duplicates
    parent_read_dedup = set(
        x for x in [*coordinates_flashed.values(), *coordinates_pe.values()]
    )
    parent_read_duplicates = parent_read_all - parent_read_dedup
    parent_read_duplicates = dict.fromkeys(parent_read_duplicates)

    with xopen.xopen(output, "w") as w:
        ujson.dump(parent_read_duplicates, w)


@reporters_deduplicate.command()
@click.argument('slices_fn')
@click.option('-d', '--duplicated_ids', help='Path to duplicated ids file')
@click.option('-o', '--output', help='Output file for deduplicated slices')
@click.option(
    "--buffer",
    help="Number of fragments to process at one time",
    default=1e6,
    type=click.INT,
)
@click.option("--stats_prefix", help="Output prefix for stats file")
@click.option("--sample_name", help="Name of sample e.g. DOX_treated_1")
@click.option("--read_type", help="Type of read", default="flashed", type=click.Choice(["flashed", "pe"], case_sensitive=False))
def remove(
    slices_fn,
    duplicated_ids,
    output="dedup.slices.tsv.gz",
    buffer=1e6,
    sample_name="",
    read_type="",
    stats_prefix="",
):

    if os.path.exists(output):
        os.unlink(output)

    df_slices = pd.read_csv(slices_fn, sep="\t", chunksize=buffer)
    ids_duplicated = set(load_json(duplicated_ids))
    n_reads_total = 0
    n_reads_unique = 0

    for ii, df in enumerate(df_slices):

        n_reads_total += df["parent_read"].nunique()

        df = (
            df.assign(parent_read_hashed=lambda df: hash_column(df["parent_read"]))
            .set_index("parent_read_hashed")
            .loc[lambda df: ~df.index.isin(ids_duplicated)]
        )

        n_reads_unique += df["parent_read"].nunique()

        df.reset_index(drop=True).to_csv(
            output, sep="\t", index=None, header=True if ii < 1 else False, mode="a"
        )

    # Sort stats
    df_stats = pd.DataFrame()
    df_stats["stat_type"] = ["not-deduplicated", "deduplicated"]
    df_stats["stat"] = [n_reads_total, n_reads_unique]
    df_stats["sample"] = sample_name
    df_stats["read_type"] = read_type
    df_stats["read_number"] = 0
    df_stats["stage"] = "deduplicate_slices"
    df_stats.to_csv(f"{stats_prefix}.read.stats.csv", index=False)