from pathlib import Path

import typer

from capcruncher.cli.common import (
    HELP_SETTINGS,
    CompressionLevelOption,
    MinimumSliceLengthOption,
    NCoresOption,
    NReadsOption,
)
from capcruncher.types import FastqSplitMethod, ReadType

fastq_app = typer.Typer(
    help="Contains methods for fastq splitting, deduplicating and digestion.",
    context_settings=HELP_SETTINGS,
    no_args_is_help=True,
)


@fastq_app.callback()
def fastq() -> None:
    """Contains methods for fastq splitting, deduplicating and digestion."""


@fastq_app.command()
def split(
    input_files: list[str] = typer.Argument(...),
    method: FastqSplitMethod = typer.Option(
        FastqSplitMethod.UNIX,
        "-m",
        "--method",
        help="Method to use for splitting.",
    ),
    output_prefix: str = typer.Option(
        "split",
        "-o",
        "--output-prefix",
        "--output_prefix",
        help="Output prefix for deduplicated fastq file(s).",
    ),
    compression_level: CompressionLevelOption = 5,
    n_reads: NReadsOption = 1_000_000,
    gzip: bool = typer.Option(
        False,
        "--gzip/--no-gzip",
        help="Determines if files are gziped or not.",
    ),
    n_cores: NCoresOption = 1,
    suffix: str = typer.Option(
        "",
        "-s",
        "--suffix",
        help="Suffix to add to output files.",
    ),
) -> None:
    """
    Splits fastq file(s) into equal chunks of n reads.
    """
    from capcruncher.api.fastq import split_fastq

    split_fastq(
        input_files=input_files,
        method=method,
        output_prefix=output_prefix,
        compression_level=compression_level,
        n_reads=n_reads,
        gzip=gzip,
        n_cores=n_cores,
        suffix=suffix,
    )


@fastq_app.command()
def digest(
    fastqs: list[str] = typer.Argument(...),
    restriction_enzyme: str = typer.Option(
        ...,
        "-r",
        "--restriction-enzyme",
        "--restriction_enzyme",
        help="Restriction enzyme name or sequence to use for in silico digestion.",
    ),
    mode: ReadType = typer.Option(
        ...,
        "-m",
        "--mode",
        help="Digestion mode. Combined (Flashed) or non-combined (PE) read pairs.",
    ),
    output_file: str = typer.Option(
        "out.fastq.gz",
        "-o",
        "--output-file",
        "--output_file",
    ),
    minimum_slice_length: MinimumSliceLengthOption = 18,
    statistics: str = typer.Option(
        "stats",
        "--statistics",
        help="Output path for stats file.",
    ),
    sample_name: str = typer.Option(
        "sampleX",
        "--sample-name",
        help="Name of sample e.g. DOX_treated_1. Required for correct statistics.",
    ),
) -> None:
    """
    Performs in silico digestion of one or a pair of fastq files.
    """
    from capcruncher.api.fastq import digest_fastq

    digest_fastq(
        fastqs=fastqs,
        restriction_site=restriction_enzyme,
        mode=mode,
        output_file=output_file,
        minimum_slice_length=minimum_slice_length,
        statistics=statistics,
        sample_name=sample_name,
    )


@fastq_app.command()
def deduplicate(
    fastq1: list[str] = typer.Option(
        ...,
        "-1",
        "--fastq1",
        help="Read 1 FASTQ files.",
    ),
    fastq2: list[str] = typer.Option(
        ...,
        "-2",
        "--fastq2",
        help="Read 2 FASTQ files.",
    ),
    output_prefix: str = typer.Option(
        "deduped",
        "-o",
        "--output-prefix",
        help="Output prefix for deduplicated FASTQ files.",
    ),
    sample_name: str = typer.Option(
        "sampleX",
        "--sample-name",
        help="Name of sample e.g. DOX_treated_1.",
    ),
    statistics: str = typer.Option(
        "stats.csv",
        "-s",
        "--statistics",
        help="Statistics output file name.",
    ),
    shuffle: bool = typer.Option(
        False,
        "--shuffle",
        help="Shuffle reads before deduplication.",
    ),
) -> None:
    """
    Identifies PCR duplicate fragments from FASTQ files.
    """
    from capcruncher.api.fastq import deduplicate_fastq

    deduplicate_fastq(
        fastq_1=[Path(fastq) for fastq in fastq1],
        fastq_2=[Path(fastq) for fastq in fastq2],
        output_prefix=output_prefix,
        statistics=statistics,
        sample_name=sample_name,
        shuffle=shuffle,
    )


cli = typer.main.get_command(fastq_app)
