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
    default="annotated.slices.parquet",
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
    "--statistics",
    help="Output path for stats file",
    default="filtering_stats.json",
)
@click.option("--sample-name", help="Name of sample e.g. DOX_treated_1")
@click.option(
    "--read-type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
@click.option(
    "--fragments/--no-fragments",
    help="Determines if read fragment aggregations are produced",
    default=True,
)
def filter(*args, **kwargs):
    """
    Removes unwanted aligned slices and identifies reporters.

    Parses a BAM file and merges this with a supplied annotation to identify unwanted slices.
    Filtering can be tuned for Capture-C, Tri-C and Tiled-C data to ensure optimal filtering.

    """
    from capcruncher.cli.alignments_filter import filter

    filter(*args, **kwargs)

