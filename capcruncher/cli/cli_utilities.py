from typing import Iterable, List, Literal
import click
from capcruncher.cli import UnsortedGroup
import ast
import pandas as pd
from cgatcore.iotools import touch_file
import os
import logging


def strip_cmdline_args(args):

    formatted_args = dict()
    for arg in args:
        key, value = arg.split("=")

        if "\\t" in value:
            formatted_args[key] = "\t"
        else:

            try:
                formatted_args[key] = ast.literal_eval(value)
            except SyntaxError:
                formatted_args[key] = value

    return formatted_args


@click.group()
def cli():
    """Contains miscellaneous functions"""


@cli.command()
@click.argument("gtf")
@click.option("-o", "--output", help="Output file name")
def gtf_to_bed12(gtf: str, output: str):
    from pybedtools import BedTool
    from capcruncher.utils import gtf_line_to_bed12_line

    bt_gtf = BedTool(gtf)
    df_gtf = bt_gtf.to_dataframe()
    df_gtf["geneid"] = df_gtf["attributes"].str.extract(r"gene_id\s?\"(.*?)\";.*")
    df_gtf = df_gtf.query('feature.isin(["5UTR", "3UTR", "exon"])')
    df_gtf = df_gtf.loc[
        lambda df: df["seqname"].str.contains(r"^chr[xXYy]?[1-9]?[0-9]?$")
    ]

    with open(output, "w") as w:
        for gene, df in df_gtf.sort_values(["seqname", "start"]).groupby("geneid"):
            w.write(gtf_line_to_bed12_line(df) + "\n")


@cli.command()
@click.argument("infiles", nargs=-1, required=True)
@click.option("-s", "--partition-size", default="2GB")
@click.option("-o", "--out-glob", required=True)
@click.option("-r", "--read-args", multiple=True)
@click.option("-w", "--write-args", multiple=True)
def repartition_csvs(
    infiles: Iterable,
    out_glob: str,
    partition_size: str,
    read_args: tuple = None,
    write_args: tuple = None,
):

    import dask.dataframe as dd

    read_args = strip_cmdline_args(read_args)
    write_args = strip_cmdline_args(write_args)

    (
        dd.read_csv(infiles, **read_args)
        .repartition(partition_size=partition_size)
        .to_csv(out_glob, **write_args)
    )


@cli.command()
@click.argument("slices")
@click.option("-o", "--output", help="Output file name")
@click.option("-m", "--method", type=click.Choice(["capture", "tri", "tiled"]))
@click.option("--file-type", type=click.Choice(["hdf5", "tsv"]))
@click.option("--sample-name", help="Name of sample e.g. DOX_treated_1")
@click.option(
    "--read-type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
def cis_and_trans_stats(
    slices: str,
    output: str,
    method: Literal["capture", "tri", "tiled"],
    file_type: str = "hdf5",
    sample_name: str = "",
    read_type: str = "",
):

    from capcruncher.tools.filter import (
        CCSliceFilter,
        TriCSliceFilter,
        TiledCSliceFilter,
    )

    filters = {
        "capture": CCSliceFilter,
        "tri": TriCSliceFilter,
        "tiled": TiledCSliceFilter,
    }
    slice_filterer = filters.get(method)

    if file_type == "tsv":
        df_slices = pd.read_csv(slices, sep="\t")

        try:
            slice_filterer(
                df_slices, sample_name=sample_name, read_type=read_type
            ).cis_or_trans_stats.to_csv(output, index=False)
        except Exception as e:
            logging.info(f"Exception: {e} occured with {slices}")
            touch_file(output)

    elif file_type == "hdf5":

        slices_iterator = pd.read_hdf(slices, "slices", chunksize=1e6, iterator=True)
        slice_stats = pd.DataFrame()

        for df in slices_iterator:
            sf = slice_filterer(df, sample_name=sample_name, read_type=read_type)
            stats = sf.cis_or_trans_stats

            slice_stats = (
                pd.concat([slice_stats, stats])
                .groupby(["viewpoint", "cis/trans", "sample", "read_type"])
                .sum()
                .reset_index()
            )
        
        slice_stats.to_csv(output, index=False)


@cli.command()
@click.argument("infiles", nargs=-1, required=True)
@click.option("-o", "--outfile", help="Output file name")
@click.option(
    "-i",
    "--index-cols",
    help="Columns to use as data_columns for queries",
    multiple=True,
)
@click.option(
    "-c",
    "--category-cols",
    help="Columns to use as data_columns for queries",
    multiple=True,
)
def merge_capcruncher_hdfs(
    infiles: Iterable,
    outfile: os.PathLike,
    index_cols: List[str] = None,
    category_cols: List[str] = None,
):

    import dask.dataframe as dd
    import dask.distributed

    client = dask.distributed.Client(n_workers=4, dashboard_address=None, processes=True)

    kwargs = {}
    kwargs.update({"data_columns": index_cols} if index_cols else {})

    ddf = dd.read_hdf(infiles, key="slices")

    # TODO: Bug when selecting columns using a string category
    # To fix, explicitly convert to string before writing
    transform_to_string = dict()
    transformed_categories = dict()
    max_len = dict()
    for col in index_cols:
        if ddf[col].dtype == "category":
            transform_to_string[col] = str
            transformed_categories[col] = list(ddf[col].cat.categories)
            max_len[col] = max([len(cat) for cat in ddf[col].cat.categories])

    ddf = ddf.astype(transform_to_string)

    if category_cols:
        ddf = ddf.categorize(columns=[*category_cols])

    # breakpoint()
    ddf = ddf.sort_values([*index_cols], npartitions="auto")

    ddf.to_hdf(
        outfile,
        key="slices",
        complib="blosc",
        complevel=2,
        min_itemsize=max_len,
        **kwargs,
    )

    with pd.HDFStore(outfile, "a") as store:
        store.create_table_index("slices", columns=index_cols, kind="full")

        if transformed_categories:
            store.put("/slices_category_metadata", pd.DataFrame(transformed_categories))

    client.shutdown()
