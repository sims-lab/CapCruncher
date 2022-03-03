from email.policy import default
import subprocess
from tempfile import NamedTemporaryFile
from typing import Iterable, List, Literal
import click
from capcruncher.cli import UnsortedGroup
import ast
import pandas as pd
from cgatcore.iotools import touch_file
import os
import logging
import glob
from capcruncher.utils import get_file_type


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
@click.option("--file-type", type=click.Choice(["parquet", "hdf5", "tsv"]))
@click.option("--sample-name", help="Name of sample e.g. DOX_treated_1")
@click.option(
    "--read-type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
@click.option(
    "-p",
    "--n_cores",
    help="Number of parallel processes to use",
    default=1,
)
@click.option(
    "--memory-limit",
    help="Maximum amount of memory to use.",
    default="1G",
)
def cis_and_trans_stats(
    slices: str,
    output: str,
    method: Literal["capture", "tri", "tiled"],
    file_type: str = "hdf5",
    sample_name: str = "",
    read_type: str = "",
    n_cores: int = 1,
    memory_limit: str = "1G",
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

    elif file_type == "parquet":

        import dask.dataframe as dd
        import dask.distributed

        with dask.distributed.Client(
            n_workers=n_cores,
            dashboard_address=None,
            processes=True,
            scheduler_port=0,
            local_directory=os.environ.get("TMPDIR", "/tmp/"),
            memory_limit=memory_limit,
        ) as client:

            ddf = dd.read_parquet(slices, engine="pyarrow")
            ddf_cis_trans_stats = ddf.map_partitions(
                lambda df: slice_filterer(
                    df, sample_name=sample_name, read_type=read_type
                ).cis_or_trans_stats
            )
            ddf_cis_trans_stats_summary = (
                ddf_cis_trans_stats.groupby(
                    ["viewpoint", "cis/trans", "sample", "read_type"]
                )
                .sum()
                .reset_index()
            )
            ddf_cis_trans_stats_summary.to_csv(output, index=False, single_file=True)


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
    help="Categorical columns",
    multiple=True,
)
@click.option(
    "-p", "--n-cores", help="Number of processes to use for merging", type=click.INT
)
def merge_capcruncher_slices(
    infiles: Iterable,
    outfile: os.PathLike,
    index_cols: List[str] = None,
    category_cols: List[str] = None,
    n_cores: int = 1,
):

    import dask.dataframe as dd
    import dask.distributed

    with dask.distributed.Client(
        n_workers=n_cores,
        dashboard_address=None,
        processes=True,
        scheduler_port=0,
        local_directory=os.environ.get("TMPDIR", "/tmp/")
    ) as client:
    
        storage_kwargs = {}
        output_format = get_file_type(outfile)

        if output_format == "hdf5":
            storage_kwargs.update({"data_columns": index_cols} if index_cols else {})

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
                **storage_kwargs,
            )

            with pd.HDFStore(outfile, "a") as store:
                store.create_table_index("slices", columns=index_cols, kind="full")

                if transformed_categories:
                    store.put(
                        "/slices_category_metadata",
                        pd.DataFrame(transformed_categories),
                    )

        elif output_format == "parquet":

            ddf = dd.read_parquet(
                infiles, chunksize="100MB", aggregate_files=True, engine="pyarrow"
            )
            ddf.to_parquet(outfile, compression="snappy", engine="pyarrow")


def dict_to_fasta(d, path):
    with open(path, "w") as fasta:
        for k, v in d.items():
            fasta.write(f">{k}\n{v}\n")

    return path


@cli.command()
@click.option("-v", "--viewpoints", help="Path to viewpoints", required=True)
@click.option("-g", "--genome", help="Path to genome fasta file", required=True)
@click.option(
    "-i", "--genome-indicies", help="Path to genome bowtie2 indices", required=True
)
@click.option("-r", "--recognition-site", help="Restriction site used", default="dpnii")
@click.option(
    "-o", "--output", help="Output file name", default="viewpoint_coordinates.bed"
)
def viewpoint_coordinates(
    viewpoints: os.PathLike,
    genome: os.PathLike,
    genome_indicies: os.PathLike = None,
    recognition_site: str = "dpnii",
    output: os.PathLike = "viewpoint_coordinates.bed",
):

    # import dask.distributed
    import concurrent.futures
    from capcruncher.cli import genome_digest
    from pybedtools import BedTool

    # client = dask.distributed.Client(
    #     n_workers=4, dashboard_address=None, processes=True
    # )

    digested_genome = NamedTemporaryFile("r+")
    viewpoints_fasta = NamedTemporaryFile("r+")
    viewpoints_aligned_bam = NamedTemporaryFile("r+")

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:

        # Digest genome to find restriction fragments
        digestion = executor.submit(
            genome_digest.digest,
            **dict(
                input_fasta=genome,
                recognition_site=recognition_site,
                output_file=digested_genome.name,
                sort=True,
            ),
        )

        # Generate a fasta file of viewpoints
        if ".fa" in viewpoints:
            fasta = viewpoints
        elif viewpoints.endswith(".tsv") or viewpoints.endswith(".csv"):
            df = pd.read_table(viewpoints)
            cols = df.columns
            fasta = dict_to_fasta(
                df.set_index(cols[0])[cols[1]].to_dict(), viewpoints_fasta.name
            )
        else:
            raise ValueError("Oligos not provided in the correct format (FASTA/TSV)")

        # Align viewpoints to the genome
        # if not genome_indicies or not os.path.exists(os.path.join(genome_indicies, ".1.bt2")):
        #     raise ValueError("No indices supplied for alignment")

        p_alignment = subprocess.Popen(
            ["bowtie2", "-x", genome_indicies, "-f", "-U", fasta],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        p_bam = subprocess.Popen(
            ["samtools", "view", "-b", "-"],
            stdout=viewpoints_aligned_bam,
            stdin=p_alignment.stdout,
        )
        p_alignment.stdout.close()
        aligned_res = p_bam.communicate()

        # Ensure genome has been digested in this time
        digestion.result()

        # Intersect digested genome with viewpoints
        bt_genome = BedTool(digested_genome.name)
        bt_viewpoints = BedTool(viewpoints_aligned_bam.name)

        intersections = bt_genome.intersect(bt_viewpoints, wa=True, wb=True)

        # Write results to file
        (
            intersections.to_dataframe()
            .drop_duplicates("name")
            .assign(oligo_name=lambda df: df["thickEnd"].str.split("_L").str[0])[
                ["chrom", "start", "end", "oligo_name"]
            ]
            .to_csv(output, index=False, header=False, sep="\t")
        )

        for tmp in [digested_genome, viewpoints_fasta, viewpoints_aligned_bam]:
            tmp.close()
