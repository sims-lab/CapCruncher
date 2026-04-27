import os
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Iterable

from loguru import logger
import typer

from capcruncher.cli.common import HELP_SETTINGS
from capcruncher.types import Assay


utilities_app = typer.Typer(
    help="Contains miscellaneous functions.",
    context_settings=HELP_SETTINGS,
    no_args_is_help=True,
)


def _first_existing_column(df: Any, candidates: Iterable[str]) -> str:
    for column in candidates:
        if column in df.columns:
            return column

    raise KeyError(
        f"None of the expected columns were present: {', '.join(candidates)}"
    )


def _has_parquet_files(path: str) -> bool:
    if os.path.isdir(path):
        return any(
            os.path.isfile(os.path.join(path, filename))
            and filename.endswith(".parquet")
            for filename in os.listdir(path)
        )

    return os.path.isfile(path)


@utilities_app.callback()
def utilities() -> None:
    """Contains miscellaneous functions"""


@utilities_app.command()
def gtf_to_bed12(
    gtf: str = typer.Argument(...),
    output: str = typer.Option(..., "-o", "--output", help="Output file name."),
) -> None:
    """
    Converts a GTF file to a BED12 file containing only 5' UTRs, 3' UTRs, and exons.

    Args:
        gtf (str): Path to the input GTF file.
        output (str): Path to the output BED12 file.

    Returns:
        None
    """

    from capcruncher.utils import gtf_line_to_bed12_line
    import pandas as pd

    gtf_cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ]
    df_gtf = pd.read_csv(gtf, sep="\t", comment="#", header=None, names=gtf_cols)
    df_gtf["geneid"] = df_gtf["attributes"].str.extract(r"gene_id\s?\"(.*?)\";.*")
    df_gtf = df_gtf.loc[df_gtf["feature"].isin(["5UTR", "3UTR", "exon"])]
    df_gtf = df_gtf.loc[df_gtf["seqname"].str.contains(r"^chr[xXYy]?[1-9]?[0-9]?$")]

    with open(output, "w") as w:
        for gene, df in df_gtf.sort_values(["seqname", "start"]).groupby("geneid"):
            w.write(gtf_line_to_bed12_line(df) + "\n")


@utilities_app.command()
def cis_and_trans_stats(
    slices: str = typer.Argument(...),
    output: str = typer.Option(..., "-o", "--output", help="Output file name."),
    sample_name: str = typer.Option(
        ..., "--sample-name", help="Name of sample e.g. DOX_treated_1."
    ),
    assay: Assay = typer.Option(
        Assay.CAPTURE, "--assay", help="Assay used to generate slices."
    ),
) -> None:
    import polars as pl
    from capcruncher.api.statistics import CisOrTransStats

    if not _has_parquet_files(slices):
        stats = CisOrTransStats(stats=[])
        with open(output, "w") as f:
            f.write(stats.model_dump_json())
        return

    parquet_path = slices
    if os.path.isdir(parquet_path):
        parquet_path = os.path.join(parquet_path, "*.parquet")

    tbl = (
        pl.scan_parquet(parquet_path)
        .with_columns(pl.col("capture").fill_null("reporter"))
        .select(["capture", "parent_id", "chrom", "viewpoint", "pe"])
    )

    if assay in {Assay.CAPTURE, Assay.TRI}:
        tbl_reporter = tbl.filter(pl.col("capture") == "reporter").select(
            ["parent_id", "chrom"]
        )
        tbl_reporter = tbl_reporter.rename({"chrom": "chrom_reporter"})

        tbl_capture = (
            tbl.filter(pl.col("capture") != "reporter")
            .select(["parent_id", "viewpoint", "chrom", "pe"])
            .rename({"chrom": "chrom_capture"})
        )

        tbl_merge = tbl_capture.join(
            tbl_reporter,
            on="parent_id",
            how="left",
        )

        tbl_merge = tbl_merge.with_columns(
            (pl.col("chrom_capture") == pl.col("chrom_reporter")).alias("is_cis")
        )

        df_cis_and_trans = tbl_merge.group_by(["viewpoint", "is_cis", "pe"]).agg(
            pl.len().alias("count")
        ).collect()

    else:
        viewpoint_chroms = (
            tbl.filter(pl.col("capture") != "reporter")
            .group_by(["viewpoint", "chrom"])
            .agg(pl.len().alias("chrom_count"))
            .sort(["viewpoint", "chrom_count"], descending=[False, True])
            .group_by("viewpoint")
            .agg(pl.col("chrom").first().alias("cis_chrom"))
        )

        df_cis_and_trans = (
            tbl.join(viewpoint_chroms, on="viewpoint", how="left")
            .with_columns((pl.col("chrom") == pl.col("cis_chrom")).alias("is_cis"))
            .group_by(["viewpoint", "parent_id", "is_cis", "pe"])
            .agg(pl.len().alias("count"))
            .group_by(["viewpoint", "is_cis", "pe"])
            .agg(pl.col("count").sum().alias("count"))
            .collect()
        )

    df_cis_and_trans = (
        df_cis_and_trans.to_pandas().rename(
            columns={"pe": "read_type", "is_cis": "cis/trans"}
        )
        .assign(
            sample=sample_name,
            **{
                "cis_or_trans": lambda df: df["cis/trans"].map(
                    {True: "cis", False: "trans"}
                )
            },
        )
        .loc[lambda df: ~df["viewpoint"].isna()]
        .sort_values(["viewpoint", "read_type", "cis_or_trans"])
    )

    stats = CisOrTransStats.from_dataframe(df_cis_and_trans)

    with open(output, "w") as f:
        f.write(stats.model_dump_json())


def dict_to_fasta(d, path):
    with open(path, "w") as fasta:
        for k, v in d.items():
            fasta.write(f">{k}\n{v}\n")

    return path


@utilities_app.command()
def viewpoint_coordinates(
    viewpoints: str = typer.Option(..., "-v", "--viewpoints", help="Path to viewpoints."),
    genome: str = typer.Option(..., "-g", "--genome", help="Path to genome fasta file."),
    genome_indicies: str = typer.Option(
        ..., "-i", "--genome-indicies", help="Path to genome bowtie2 indices."
    ),
    recognition_site: str = typer.Option(
        "dpnii", "-r", "--recognition-site", help="Restriction site used."
    ),
    output: str = typer.Option(
        "viewpoint_coordinates.bed", "-o", "--output", help="Output file name."
    ),
) -> None:
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

    from capcruncher.api.genome import digest_genome
    import pandas as pd
    import pyranges1 as pr
    import pysam

    def bam_to_bed_df(bam_path: os.PathLike):
        rows = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.reference_name is None:
                    continue

                rows.append(
                    {
                        "Chromosome": read.reference_name,
                        "Start": read.reference_start,
                        "End": read.reference_end,
                        "Name": read.query_name,
                    }
                )

        return pd.DataFrame.from_records(
            rows, columns=["Chromosome", "Start", "End", "Name"]
        )

    digested_genome = NamedTemporaryFile("r+")
    viewpoints_fasta = NamedTemporaryFile("r+")
    viewpoints_aligned_bam = NamedTemporaryFile("r+")

    digest_genome(
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
    if p_alignment.stdout is not None:
        p_alignment.stdout.close()
    aligned_res = p_bam.communicate()

    # Intersect digested genome with viewpoints
    gr_genome = pr.PyRanges(
        pd.read_csv(
            digested_genome.name,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Name"],
        )
    )
    gr_viewpoints = pr.PyRanges(bam_to_bed_df(Path(viewpoints_aligned_bam.name)))
    intersections = gr_genome.join_overlaps(
        gr_viewpoints, suffix="_vp", strand_behavior="ignore"
    )
    fragment_name_column = _first_existing_column(intersections, ["Name"])
    viewpoint_name_column = _first_existing_column(
        intersections, ["Name_vp", "Name_b"]
    )

    # Write results to file
    (
        intersections.drop_duplicates(fragment_name_column)
        .assign(
            oligo_name=lambda df: df[viewpoint_name_column]
            .astype(str)
            .str.split("_L")
            .str[0]
        )[["Chromosome", "Start", "End", "oligo_name"]]
        .rename(
            columns={
                "Chromosome": "chrom",
                "Start": "start",
                "End": "end",
            }
        )
        .to_csv(output, index=False, header=False, sep="\t")
    )

    for tmp in [digested_genome, viewpoints_fasta, viewpoints_aligned_bam]:
        tmp.close()


def dump_cooler(path: str, viewpoint: str, resolution: int | None = None):
    import cooler.api as cooler

    if viewpoint is None:
        raise ValueError("A viewpoint is required when dumping cooler files.")

    if resolution:
        path = cooler.Cooler(f"{path}::{viewpoint}/resolutions/{resolution}")
    else:
        path = cooler.Cooler(f"{path}::{viewpoint}")

    pixels = path.pixels()[:]
    return pixels


def dump_capcruncher_parquet(path: str, viewpoint: str | None = None):
    import polars as pl

    parquet_path = path
    if os.path.isdir(parquet_path):
        parquet_path = os.path.join(parquet_path, "*.parquet")

    if viewpoint:
        tbl = pl.scan_parquet(parquet_path).filter(pl.col("viewpoint") == viewpoint)
    else:
        tbl = pl.scan_parquet(parquet_path)

    return tbl.collect().to_pandas()


@utilities_app.command()
def dump(
    path: str = typer.Argument(...),
    viewpoint: str | None = typer.Option(
        None, "-v", "--viewpoint", help="Viewpoint to extract."
    ),
    resolution: int | None = typer.Option(
        None,
        "-r",
        "--resolution",
        help="Resolution to extract. Only used for cooler (hdf5) files.",
    ),
    output: str = typer.Option(
        "capcruncher_dump.tsv", "-o", "--output", help="Output file name."
    ),
) -> None:
    """
    Dumps the contents of a cooler or capcruncher parquet file to a TSV file

    Args:
        path (str): Path to cooler or capcruncher parquet file
        viewpoint (str, optional): Viewpoint to extract. Defaults to None.
        resolution (int, optional): Resolution to extract. Only used for cooler (hdf5) files. Defaults to None.
        output (str, optional): Output file name. Defaults to "capcruncher_dump.tsv".
    """

    assert os.path.exists(path), "File does not exist"

    if ".hdf5" in path:
        if viewpoint is None:
            raise typer.BadParameter("A viewpoint is required for cooler files.")
        df = dump_cooler(path, viewpoint, resolution)
    elif ".parquet" in path:
        df = dump_capcruncher_parquet(path, viewpoint)
    else:
        raise ValueError("File type not supported")

    df.to_csv(output, sep="\t", index=False)


@utilities_app.command()
def regenerate_fastq(
    fastq1: str = typer.Option(..., "-1", "--fastq1", help="Path to FASTQ file 1."),
    fastq2: str = typer.Option(..., "-2", "--fastq2", help="Path to FASTQ file 2."),
    parquet_file: str = typer.Option(
        ...,
        "-p",
        "--parquet-file",
        help="Path to parquet file from which to extract the required reads.",
    ),
    output_prefix: str = typer.Option(
        "regenerated_", "-o", "--output-prefix", help="Output file prefix."
    ),
) -> None:
    """
    Regenerates a FASTQ file from a parquet file containing the required reads

    Args:
        fastq1 (str): Path to the first FASTQ file
        fastq2 (str): Path to the second FASTQ file
        parquet_file (str, optional): Path to the parquet file from which to extract the required reads. Defaults to None.
        output (str, optional): Prefix for the output file. Defaults to "regenerated_".

    Raises:
        AssertionError: If the specified parquet file does not exist.

    Returns:
        None
    """
    import pathlib

    import polars as pl
    import pysam
    from xopen import xopen

    assert os.path.exists(parquet_file), f"File {parquet_file} does not exist"

    parquet_file_path = pathlib.Path(parquet_file)
    if parquet_file_path.is_dir():
        parquet_file = str(parquet_file_path / "*.parquet")

    outpath = pathlib.Path(output_prefix).with_suffix("")

    logger.info(f"Extracting reads info from {parquet_file}")
    with pl.StringCache():
        read_names = set(
            pl.scan_parquet(parquet_file)
            .select("parent_read")
            .unique()
            .collect()["parent_read"]
            .to_list()
        )

        logger.info(f"Writing reads to {outpath}")
        with pysam.FastxFile(fastq1) as r1:
            with pysam.FastxFile(fastq2) as r2:
                with xopen(f"{outpath}_1.fastq.gz", "w") as w1:
                    with xopen(f"{outpath}_2.fastq.gz", "w") as w2:
                        for read_1, read_2 in zip(r1, r2):
                            if read_1.name in read_names:
                                w1.write(str(read_1) + "\n")
                                w2.write(str(read_2) + "\n")

    logger.info("Done")


@utilities_app.command()
def make_chicago_maps(
    fragments: str = typer.Option(
        "capcruncher_output/resources/restriction_fragments/genome.digest.bed.gz",
        "--fragments",
        help="Path to fragments file.",
    ),
    viewpoints: str = typer.Option(
        ..., "--viewpoints", help="Path to viewpoints file used for capcruncher."
    ),
    outputdir: str = typer.Option(..., "-o", "--outputdir", help="Path to output directory."),
) -> None:
    """
    Restriction map file (.rmap) - a bed file containing coordinates of the restriction fragments. By default, 4 columns: chr, start, end, fragmentID.
    Bait map file (.baitmap) - a bed file containing coordinates of the baited restriction fragments, and their associated annotations. By default, 5 columns: chr, start, end, fragmentID, baitAnnotation. The regions specified in this file, including their fragmentIDs, must be an exact subset of those in the .rmap file. The baitAnnotation is a text field that is used only to annotate the output and plots.
    """
    import pathlib
    import pyranges1 as pr

    # Rename fragments file to suit chicago
    fragments_new = pathlib.Path(outputdir) / (pathlib.Path(fragments).stem + ".rmap")
    if not fragments_new.exists():
        fragments_new.symlink_to(pathlib.Path(fragments).resolve())

    # Baitmap file
    viewpoints = pr.read_bed(viewpoints)
    fragments = pr.read_bed(fragments)

    intersections = fragments.join_overlaps(
        viewpoints, suffix="_vp", strand_behavior="ignore"
    )
    fragment_name_column = _first_existing_column(intersections, ["Name"])
    viewpoint_name_column = _first_existing_column(
        intersections, ["Name_vp", "Name_b"]
    )

    df_baitmap = (
        intersections[
            ["Chromosome", "Start", "End", fragment_name_column, viewpoint_name_column]
        ]
        .rename(
            columns={
                "Chromosome": "chr",
                "Start": "start",
                "End": "end",
                fragment_name_column: "baitAnnotation",
                viewpoint_name_column: "fragmentID",
            }
        )
    )

    df_baitmap.to_csv(
        os.path.join(outputdir, "viewpoints.baitmap"),
        sep="\t",
        index=False,
        header=False,
    )


cli = typer.main.get_command(utilities_app)
