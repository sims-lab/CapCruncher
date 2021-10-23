from typing import Iterable, Literal
import click
from capcruncher.cli import UnsortedGroup
import ast


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
@click.option("--sample_name", help="Name of sample e.g. DOX_treated_1")
@click.option(
    "--read_type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
def get_cis_and_trans_stats(
    slices: str,
    output: str,
    method: Literal["capture", "tri", "tiled"],
    sample_name: str = "",
    read_type: str = "",
):

    import pandas as pd
    from cgatcore.iotools import touch_file
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

    df_slices = pd.read_csv(slices, sep="\t")

    try:
        slice_filterer(
            df_slices, sample_name=sample_name, read_type=read_type
        ).cis_or_trans_stats.to_csv(output, index=False)
    except:
        touch_file(output)
