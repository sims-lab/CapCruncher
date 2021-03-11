import pandas as pd

import click
import xopen
from ccanalyser.cli._reporters import cli
from ccanalyser.utils import hash_column, load_json
import ujson
import os


@cli.group()
def deduplicate():
    """Identifies duplicate aligned fragments and removes them."""


@deduplicate.command()
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
@click.option(
    "--read_type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
def identify(
    fragments_fn: os.PathLike,
    output: os.PathLike = "duplicated_ids.json",
    buffer: int = 1e6,
    read_type: str = "flashed",
):
    """Identifies read fragments with duplicate coordinates and outputs these to a .json file.

    Args:
     fragments_fn (os.PathLike): Input fragments.tsv file to process.
     output (os.PathLike, optional): Output path to output duplicated parental read ids. Defaults to "duplicated_ids.json".
     buffer (int, optional): Number of fragments to process in memory. Defaults to 1e6.
     read_type (str, optional): Process combined(flashed) or non-combined reads (pe). 
                                Due to the low confidence in the quaility of pe reads, duplicates are identified by
                                removing any fragments with matching start and end coordinates.
                                Defaults to "flashed".
    """

    fragments = pd.read_csv(
        fragments_fn,
        sep="\t",
        chunksize=buffer,
        usecols=["parent_read", "coordinates"],
    )

    coordinates_deduplicated = dict() #{coord hash: id hash}
    fragments_all = set() # {id hash}

    for df in fragments:

        df = df.sample(
            frac=1
        )  # Shuffles to stop fragments at the end of the sample always being removed

        if read_type == "flashed":

            coords_hashed = hash_column(df["coordinates"])
            parent_read_hashed = hash_column(df["parent_read"])

        else:
            #TODO: Slight bug here as flashed coordinates will not be compared to pe.

            # Extract chrom1 + start1 + chrom(last entry) + end(last entry)
            coords_df = df["coordinates"].str.extract(
                r"^chr(?P<chrom1>[\d|X|Y|M]+):(?P<start>\d+).*\|chr(?P<chrom2>[\d|X|Y|M]+):\d+-(?P<end>\d+)"
            )

            # {chrom1+start1+chrom-1+end-1(hashed): id(hashed)}
            coords_hashed = hash_column(
                coords_df["chrom1"].str.cat(coords_df.iloc[:, 1:])
            )
            parent_read_hashed = hash_column(df["parent_read"])

        coordinates_deduplicated_sample = dict(zip(coords_hashed, parent_read_hashed))
        coordinates_deduplicated.update(coordinates_deduplicated_sample) # Use dict to remove duplicated coords
        fragments_all.update({x for x in parent_read_hashed}) # Store all ids for duplicate identification

    # Identify duplicates
    fragments_no_dup = {x for x in coordinates_deduplicated.values()}
    fragments_dup = fragments_all - fragments_no_dup

    with xopen.xopen(output, "w") as w:
        ujson.dump(dict.fromkeys(fragments_dup), w)


@deduplicate.command()
@click.argument("slices_fn")
@click.option("-d", "--duplicated_ids", help="Path to duplicated ids file")
@click.option("-o", "--output", help="Output file for deduplicated slices")
@click.option(
    "--buffer",
    help="Number of fragments to process at one time",
    default=1e6,
    type=click.INT,
)
@click.option("--stats_prefix", help="Output prefix for stats file")
@click.option("--sample_name", help="Name of sample e.g. DOX_treated_1")
@click.option(
    "--read_type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
def remove(
    slices_fn: os.PathLike,
    duplicated_ids: os.PathLike,
    output: os.PathLike = "dedup.slices.tsv.gz",
    buffer: int = 1e6,
    sample_name: str = "",
    read_type: str = "",
    stats_prefix: os.PathLike = "",
):
    """Removes duplicated fragments from slices file according to duplicated reads found by identify.


    Args:
     slices_fn (os.PathLike): Input slices.tsv file.
     duplicated_ids (os.PathLike): Duplicated parental read ids in json format.
     output (os.PathLike, optional): Output file path for deduplicated slices. Defaults to "dedup.slices.tsv.gz".
     buffer (int, optional): Number of slices to process in memory. Defaults to 1e6.
     sample_name (str, optional): Name of sample being processed e.g. DOX-treated_1 used for statistics. Defaults to "".
     read_type (str, optional): Process combined(flashed) or non-combined reads (pe) used for statistics. Defaults to "".
     stats_prefix (os.PathLike, optional): Output path for deduplication statistics. Defaults to "".
    """

    if os.path.exists(output):
        os.unlink(output)

    df_slices = pd.read_csv(slices_fn, sep="\t", chunksize=buffer)
    ids_duplicated = {int(x) for x in load_json(duplicated_ids)}
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