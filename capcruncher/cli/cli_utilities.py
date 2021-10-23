from typing import Iterable
import click
from capcruncher.cli import UnsortedGroup


def strip_cmdline_args(args):

    formatted_args = dict()
    try:
        for arg in args.split():
            if arg:
                key, value = arg.split("=")
                formatted_args[key] = value if not value == "None" else None
    except AttributeError as e:
        pass

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
@click.option("--read-args")
@click.option("--write-args")
def repartition_csvs(
    infiles: Iterable,
    out_glob: str,
    partition_size: str,
    read_args=None,
    write_args=None,
):

    import dask.dataframe as dd

    read_args = strip_cmdline_args(read_args)
    write_args = strip_cmdline_args(write_args)

    (
        dd.read_csv(infiles, **read_args)
        .repartition(partition_size=partition_size)
        .to_csv(out_glob, **write_args)
    )
