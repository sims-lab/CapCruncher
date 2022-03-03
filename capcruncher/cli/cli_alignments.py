import click


@click.group()
def cli():
    """Contains methods for reporter annotating, identifying and deduplication."""


@cli.command()
@click.argument("slices")
@click.option(
    "-a",
    "--actions",
    help="Determines if the overlaps are counted or if the name should just be reported",
    multiple=True,
    type=click.Choice(
        ["get", "count"],
    ),
)
@click.option(
    "-b", "--bed_files", help="Bed file(s) to intersect with slices", multiple=True
)
@click.option(
    "-n",
    "--names",
    help="Names to use as column names for the output tsv file.",
    multiple=True,
)
@click.option(
    "-f",
    "--overlap_fractions",
    help="The minimum overlap required for an intersection between two intervals to be reported.",
    multiple=True,
    default=[
        1e-9,
    ],
    type=click.FLOAT,
)
@click.option(
    "-t",
    "--dtypes",
    help="Data type for column",
    multiple=True,
    default=[
        "str",
    ],
)
@click.option(
    "-o",
    "--output",
    help="Path for the annotated slices to be output.",
    default="annotated.slices.tsv.gz",
)
@click.option(
    "--duplicates",
    help="Method to use for reconciling duplicate slices (i.e. multimapping). Currently only 'remove' is supported.",
    type=click.Choice(["remove"]),
    default="remove",
)
@click.option(
    "-p",
    "--n_cores",
    help="Intersections are performed in parallel, set this to the number of intersections required",
    default=1,
)
@click.option(
    "--invalid_bed_action",
    help=" ".join(
        [
            "Method to deal with invalid bed files e.g. blank or incorrectly formatted.",
            "Setting this to 'ignore' will report default N/A values (either '.' or 0) for invalid files",
        ]
    ),
    default="error",
    type=click.Choice(["ignore", "error"]),
)
@click.option(
    "--blacklist",
    help="Regions to remove from the BAM file prior to annotation",
)
@click.option(
    "--prioritize-cis-slices",
    is_flag=True,
    help="Attempts to prevent slices on the most common chromosome in a fragment (ideally cis to the viewpoint) being removed by deduplication",
)
@click.option(
    "--priority-chroms",
    help="A comma separated list of chromosomes to prioritize during deduplication",
)
def annotate(*args, **kwargs):
    """
    Annotates a bed file with other bed files using bedtools intersect.

    Whilst bedtools intersect allows for interval names and counts to be used for annotating intervals, this command
    provides the ability to annotate intervals with both interval names and counts at the same time. As the pipeline allows
    for empty bed files, this command has built in support to deal with blank/malformed bed files and will return default N/A values.

    Prior to interval annotation, the bed file to be intersected is validated and duplicate entries/multimapping reads are removed
    to ensure consistent annotations and prevent issues with reporter identification.

    """

    from capcruncher.cli.alignments_annotate import annotate

    annotate(*args, **kwargs)


@cli.command()
@click.argument("method", type=click.Choice(["capture", "tri", "tiled"]))
@click.option("-b", "--bam", help="Bam file to process", required=True)
@click.option(
    "-a",
    "--annotations",
    help="Annotations for the bam file that must contain the required columns, see description.",
    required=True,
)
@click.option(
    "--custom-filtering",
    help="Custom filtering to be used. This must be supplied as a path to a yaml file.",
    default=None,
)
@click.option(
    "-o",
    "--output_prefix",
    help="Output prefix for deduplicated fastq file(s)",
    default="",
)
@click.option(
    "--stats-prefix",
    help="Output prefix for stats file(s)",
    default="",
)
@click.option("--sample-name", help="Name of sample e.g. DOX_treated_1")
@click.option(
    "--read-type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
# @click.option(
#     "--gzip/--no-gzip", help="Determines if files are gziped or not", default=False
# )
@click.option(
    "--fragments/--no-fragments",
    help="Determines if read fragment aggregations are produced",
    default=True,
)
@click.option(
    "--read-stats/--no-read-stats",
    help="Determines if read level statistics are output",
    default=True,
)
@click.option(
    "--slice-stats/--no-slice-stats",
    help="Determines if slice level statistics are output",
    default=True,
)
@click.option(
    "--cis-and-trans-stats/--no-cis-and-trans-stats",
    help="Determines cis/trans statistics are output",
    default=True,
)
@click.option(
    "--output-format",
    help="Determines file output format",
    default="parquet",
    type=click.Choice(["tsv", "hdf5", "parquet"]),
)
def filter(*args, **kwargs):
    """
    Removes unwanted aligned slices and identifies reporters.

    Parses a BAM file and merges this with a supplied annotation to identify unwanted slices.
    Filtering can be tuned for Capture-C, Tri-C and Tiled-C data to ensure optimal filtering.

    """
    from capcruncher.cli.alignments_filter import filter

    filter(*args, **kwargs)


@cli.group()
def deduplicate():
    """
    Identifies and removes duplicated aligned fragments.

    PCR duplicates are very commonly present in Capture-C/Tri-C/Tiled-C data and must be removed
    for accurate analysis. Unlike fastq deduplicate, this command removes fragments with identical
    genomic coordinates.

    Non-combined (pe) and combined (flashed) reads are treated slightly differently due to the increased
    confidence that the ligation junction has been captured for the flashed reads.

    """


@deduplicate.command()
@click.argument("fragments", nargs=-1, required=True)
@click.option(
    "-v",
    "--viewpoint",
    help="Viewpoint to process, leave blank for all viewpoints",
    default="",
)
@click.option(
    "-t",
    "--file-type",
    help="File format for input",
    default="auto",
    type=click.Choice(["auto", "tsv", "hdf5", "parquet"], case_sensitive=False),
)
@click.option(
    "-o",
    "--output",
    help="Path for outputting fragments with duplicated coordinates in json format.",
    default="duplicated_ids.pickle",
)
@click.option(
    "--buffer",
    help="Number of fragments to process at one time in order to preserve memory.",
    default=1e6,
    type=click.INT,
)
@click.option(
    "--read-type",
    help="Indicates if the fragments have been combined (flashed) or not (pe).",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
@click.option(
    "-p",
    "--n_cores",
    help="Number of parallel processes to use for deduplication",
    default=1,
)
@click.option(
    "--memory-limit",
    help="Maximum amount of memory to use.",
    default="1G",
)
def identify(*args, **kwargs):
    from capcruncher.cli.alignments_deduplicate import identify

    identify(*args, **kwargs)


@deduplicate.command()
@click.argument("slices", nargs=-1, required=True)
@click.option(
    "-d",
    "--duplicated_ids",
    help="Path to duplicated fragment ids determined by the 'identify' subcommand.",
)
@click.option(
    "-o",
    "--output",
    help="Path for outputting deduplicated slices in tsv format.",
    default="slices_dedup.hdf5",
)
@click.option(
    "-t",
    "--file-type",
    help="File format for input",
    default="auto",
    type=click.Choice(["auto", "tsv", "hdf5", "parquet"], case_sensitive=False),
)
@click.option(
    "--buffer",
    help="Number of fragments to process at one time, in order to preserve memory.",
    default=1e6,
    type=click.INT,
)
@click.option("--stats-prefix", help="Output prefix for deduplication statistics")
@click.option(
    "--sample-name",
    help="Name of sample being analysed e.g. DOX_treated_1. Required for correct statistics.",
)
@click.option(
    "--read-type",
    help="Indicates if the fragments have been combined (flashed) or not (pe). Required for correct statistics.",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
@click.option(
    "-p",
    "--n_cores",
    help="Number of parallel processes to use for deduplication",
    default=1,
)
@click.option(
    "--memory-limit",
    help="Maximum amount of memory to use.",
    default="1G",
)
def remove(*args, **kwargs):
    """
    Removes duplicated aligned fragments.

    Parses a tsv file containing aligned read slices and outputs only slices from unique fragments.
    Duplicated parental read id determined by the "identify" subcommand are located within the
    slices tsv file and removed.

    Outputs statistics for the number of unique slices and the number of duplicate slices identified.
    """
    from capcruncher.cli.alignments_deduplicate import remove

    remove(*args, **kwargs)
