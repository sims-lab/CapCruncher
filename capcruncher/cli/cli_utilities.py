import subprocess
from tempfile import NamedTemporaryFile
from typing import Iterable, List, Literal
import click
import pandas as pd
import os
from loguru import logger
from capcruncher.utils import get_file_type
import ibis
from ibis import _


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
@click.argument("slices")
@click.option("-o", "--output", help="Output file name")
@click.option("--sample-name", help="Name of sample e.g. DOX_treated_1")
@click.option(
    "--assay",
    help="Assay used to generate slices",
    type=click.Choice(["capture", "tri", "tiled"]),
)
def cis_and_trans_stats(
    slices: str,
    output: str,
    sample_name: str,
    assay: Literal["capture", "tri", "tiled"] = "capture",
):
    con = ibis.duckdb.connect()

    if not os.path.isdir(slices):
        tbl = con.register(f"parquet://{slices}", table_name="slices_tbl")
    else:
        tbl = con.register(f"parquet://{slices}/*.parquet", table_name="slices_tbl")

    tbl = tbl.mutate(capture=tbl["capture"].fillna("reporter")).select(
        ["capture", "parent_id", "chrom", "viewpoint", "pe"]
    )

    if assay in ["capture", "tri"]:
        tbl_reporter = tbl[(tbl["capture"] == "reporter")].drop(
            "viewpoint", "pe", "capture"
        )

        tbl_capture = tbl[~(tbl["capture"] == "reporter")]

        tbl_merge = tbl_capture.join(
            tbl_reporter,
            predicates=[
                "parent_id",
            ],
            lname="{name}_capture",
            rname="{name}_reporter",
            how="left",
        )

        tbl_merge = tbl_merge.mutate(
            is_cis=(tbl_merge["chrom_capture"] == tbl_merge["chrom_reporter"])
        )

        df_cis_and_trans = (
            tbl_merge.group_by(["viewpoint", "is_cis", "pe"])
            .aggregate(
                count=_.count(),
            )
            .execute(limit=None)
        )

    else:
        viewpoint_chroms = (
            tbl.filter(tbl["capture"] != "reporter")
            .group_by(["viewpoint", "chrom"])
            .aggregate(chrom_count=_.count())
            .order_by(ibis.desc("chrom_count"))
            .to_pandas()
            .drop_duplicates("viewpoint")
            .set_index("viewpoint")
            .to_dict()["chrom"]
        )

        chrom_mapping_exp = ibis.case()
        for k, v in viewpoint_chroms.items():
            chrom_mapping_exp = chrom_mapping_exp.when(tbl.viewpoint == k, v)
        chrom_mapping_exp = chrom_mapping_exp.end()

        df_cis_and_trans = (
            tbl.mutate(cis_chrom=chrom_mapping_exp)
            .mutate(is_cis=_.chrom == _.cis_chrom)
            .group_by(["viewpoint", "parent_id", "is_cis", "pe"])
            .aggregate(count=_.count())
            .group_by(["viewpoint", "is_cis", "pe"])
            .aggregate(count=_["count"].sum())
            .execute(limit=None)
        )

    df_cis_and_trans = (
        df_cis_and_trans.rename(columns={"pe": "read_type", "is_cis": "cis/trans"})
        .assign(
            sample=sample_name,
            **{
                "cis/trans": lambda df: df["cis/trans"].map(
                    {True: "cis", False: "trans"}
                )
            },
        )
        .loc[lambda df: ~df["viewpoint"].isna()]
        .sort_values(["viewpoint", "read_type", "cis/trans"])
    )

    df_cis_and_trans.to_csv(output, index=False)


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
        local_directory=os.environ.get("TMPDIR", "/tmp/"),
    ) as _client:
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
            import pyarrow.dataset as ds

            datasets = ds.dataset([ds.dataset(fn) for fn in infiles])
            sc = datasets.scanner()
            ds.write_dataset(
                sc,
                outfile,
                format="parquet",
                partitioning_flavor="hive",
                max_rows_per_file=1024 * 1024 * 10,
            )


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
    """
    Aligns viewpoints to a genome and returns the coordinates of the viewpoint
    in the genome.

    Viewpoints can be supplied as a FASTA file or a TSV file with the first column
    containing the name of the viewpoint and the second column containing the
    sequence of the viewpoint.

    Args:
        viewpoints (os.PathLike): Path to viewpoints
        genome (os.PathLike): Path to genome fasta file
        genome_indicies (os.PathLike, optional): Path to genome bowtie2 indices. Defaults to None.
        recognition_site (str, optional): Restriction site used. Defaults to "dpnii".
        output (os.PathLike, optional): Output file name. Defaults to "viewpoint_coordinates.bed".

    Raises:
        ValueError: If viewpoints are not supplied in the correct format
        ValueError: If no bowtie2 indices are supplied
    """

    from capcruncher.cli import genome_digest
    from pybedtools import BedTool

    digested_genome = NamedTemporaryFile("r+")
    viewpoints_fasta = NamedTemporaryFile("r+")
    viewpoints_aligned_bam = NamedTemporaryFile("r+")

    genome_digest.digest(
        input_fasta=genome,
        recognition_site=recognition_site,
        output_file=digested_genome.name,
        sort=True,
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


def dump_cooler(path: str, viewpoint: str, resolution: int = None) -> pd.DataFrame:
    import cooler.api as cooler

    if resolution:
        path = cooler.Cooler(f"{path}::{viewpoint}/resolutions/{resolution}")
    else:
        path = cooler.Cooler(f"{path}::{viewpoint}")

    pixels = path.pixels()[:]
    return pixels


def dump_capcruncher_parquet(path: str, viewpoint: str = None) -> pd.DataFrame:
    import ibis

    con = ibis.duckdb.connect()

    if viewpoint:
        tbl = con.register(f"parquet://{path}/*.parquet", table_name="slices_tbl")
        tbl = tbl.filter(tbl.viewpoint == viewpoint)
    else:
        tbl = con.register(f"parquet://{path}", table_name="slices_tbl")

    return tbl.execute(limit=None)


@cli.command()
@click.argument("path")
@click.option("-v", "--viewpoint", help="Viewpoint to extract")
@click.option(
    "-r",
    "--resolution",
    help="Resolution to extract. Only used for cooler (hdf5) files",
)
@click.option("-o", "--output", help="Output file name", default="capcruncher_dump.tsv")
def dump(
    path: str,
    viewpoint: str = None,
    resolution: int = None,
    output: str = "capcruncher_dump.tsv",
):
    """
    Dumps the contents of a cooler or capcruncher parquet file to a TSV file

    Args:
        path (str): Path to cooler or capcruncher parquet file
        viewpoint (str, optional): Viewpoint to extract. Defaults to None.
        resolution (int, optional): Resolution to extract. Only used for cooler (hdf5) files. Defaults to None.
        output (str, optional): Output file name. Defaults to "capcruncher_dump.tsv".
    """

    import pandas as pd

    assert os.path.exists(path), "File does not exist"

    if ".hdf5" in path:
        df = dump_cooler(path, viewpoint, resolution)
    elif ".parquet" in path:
        df = dump_capcruncher_parquet(path, viewpoint)
    else:
        raise ValueError("File type not supported")

    df.to_csv(output, sep="\t", index=False)
