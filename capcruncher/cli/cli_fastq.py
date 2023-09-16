import click
from capcruncher.cli import UnsortedGroup


@click.group()
def cli():
    """Contains methods for fastq splitting, deduplicating and digestion."""


@cli.command()
@click.argument("input_files", nargs=-1, required=True)
@click.option(
    "-m",
    "--method",
    help="Method to use for splitting",
    type=click.Choice(["python", "unix"]),
    default="unix",
)
@click.option(
    "-o",
    "--output_prefix",
    help="Output prefix for deduplicated fastq file(s)",
    default="split",
)
@click.option(
    "--compression_level",
    help="Level of compression for output files",
    default=5,
    type=click.INT,
)
@click.option(
    "-n",
    "--n_reads",
    help="Number of reads per fastq file",
    default=1e6,
    type=click.INT,
)
@click.option(
    "--gzip/--no-gzip", help="Determines if files are gziped or not", default=False
)
@click.option("-p", "--n_cores", default=1, type=click.INT)
@click.option(
    "-s",
    "--suffix",
    help="Suffix to add to output files (ignore {read_number}.fastq as this is added automatically)",
    default="",
)
def split(*args, **kwargs):
    """
    Splits fastq file(s) into equal chunks of n reads.

    """

    from capcruncher.cli.fastq_split import split

    split(*args, **kwargs)


@cli.command()
@click.argument("input_fastq", nargs=-1, required=True)
@click.option(
    "-r",
    "--restriction_enzyme",
    help="Restriction enzyme name or sequence to use for in silico digestion.",
    required=True,
)
@click.option(
    "-m",
    "--mode",
    help="Digestion mode. Combined (Flashed) or non-combined (PE) read pairs.",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
    required=True,
)
@click.option("-o", "--output_file", default="out.fastq.gz")
@click.option("-p", "--n_cores", default=1, type=click.INT)
@click.option("--minimum_slice_length", default=18, type=click.INT)
@click.option("--keep_cutsite", default=False)
@click.option(
    "--compression_level",
    help="Level of compression for output files (1=low, 9=high)",
    default=5,
    type=click.INT,
)
@click.option(
    "--read_buffer",
    help="Number of reads to process before writing to file to conserve memory.",
    default=1e5,
    type=click.INT,
)
@click.option("--stats-prefix", help="Output prefix for stats file", default="stats")
@click.option(
    "--sample-name",
    help="Name of sample e.g. DOX_treated_1. Required for correct statistics.",
    default="sampleX",
)
def digest(*args, **kwargs):
    """
    Performs in silico digestion of one or a pair of fastq files.
    """
    from capcruncher.cli.fastq_digest import digest

    digest(*args, **kwargs)


@cli.command()
def deduplicate():
    """
    Identifies PCR duplicate fragments from Fastq files.

    PCR duplicates are very commonly present in Capture-C/Tri-C/Tiled-C data and must be removed
    for accurate analysis. These commands attempt to identify and remove duplicate reads/fragments
    from fastq file(s) to speed up downstream analysis.

    """
