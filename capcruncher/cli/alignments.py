import typer

from capcruncher.cli.common import HELP_SETTINGS, NCoresOption


alignments_app = typer.Typer(
    help="Contains methods for reporter annotating, identifying and deduplication.",
    context_settings=HELP_SETTINGS,
    no_args_is_help=True,
)


@alignments_app.callback()
def alignments():
    """Contains methods for reporter annotating, identifying and deduplication."""


@alignments_app.command()
def annotate(
    slices: str = typer.Argument(...),
    actions: list[str] | None = typer.Option(
        None,
        "-a",
        "--actions",
        help="Determines if the overlaps are counted or if the name should just be reported.",
    ),
    bed_files: list[str] | None = typer.Option(
        None,
        "-b",
        "--bed-files",
        "--bed_files",
        help="Bed file(s) to intersect with slices.",
    ),
    names: list[str] | None = typer.Option(
        None,
        "-n",
        "--names",
        help="Names to use as column names for the output tsv file.",
    ),
    overlap_fractions: list[float] | None = typer.Option(
        None,
        "-f",
        "--overlap-fractions",
        "--overlap_fractions",
        help="The minimum overlap required for an intersection between two intervals to be reported.",
    ),
    dtypes: list[str] | None = typer.Option(
        None,
        "-t",
        "--dtypes",
        help="Data type for column.",
    ),
    output: str = typer.Option(
        "annotated.slices.parquet",
        "-o",
        "--output",
        help="Path for the annotated slices to be output.",
    ),
    duplicates: str = typer.Option(
        "remove",
        "--duplicates",
        help="Method to use for reconciling duplicate slices.",
    ),
    n_cores: NCoresOption = 1,
    invalid_bed_action: str = typer.Option(
        "error",
        "--invalid-bed-action",
        "--invalid_bed_action",
        help="Method to deal with invalid bed files.",
    ),
    blacklist: str = typer.Option(
        "",
        "--blacklist",
        help="Regions to remove from the BAM file prior to annotation.",
    ),
    prioritize_cis_slices: bool = typer.Option(
        False,
        "--prioritize-cis-slices",
        help="Attempts to prevent cis slices being removed by deduplication.",
    ),
    priority_chroms: str = typer.Option(
        "",
        "--priority-chroms",
        help="A comma separated list of chromosomes to prioritize during deduplication.",
    ),
):
    """Annotate a bed file with other bed files."""

    from capcruncher.api.alignments_annotate import annotate as annotate_alignments

    annotate_alignments(
        slices=slices,
        actions=tuple(actions or ()),
        bed_files=tuple(bed_files or ()),
        names=tuple(names or ()),
        overlap_fractions=tuple(overlap_fractions or (1e-9,)),
        dtypes=tuple(dtypes or ("str",)),
        output=output,
        duplicates=duplicates,
        n_cores=n_cores,
        invalid_bed_action=invalid_bed_action,
        blacklist=blacklist,
        prioritize_cis_slices=prioritize_cis_slices,
        priority_chroms=priority_chroms,
    )


@alignments_app.command("filter")
def filter_alignments(
    method: str = typer.Argument(..., help="Filtering method: capture, tri, or tiled."),
    bam: str = typer.Option(
        ...,
        "-b",
        "--bam",
        help="Bam file to process.",
    ),
    annotations: str = typer.Option(
        ...,
        "-a",
        "--annotations",
        help="Annotations for the bam file.",
    ),
    custom_filtering: str | None = typer.Option(
        None,
        "--custom-filtering",
        help="Custom filtering yaml file.",
    ),
    output_prefix: str = typer.Option(
        "",
        "-o",
        "--output-prefix",
        "--output_prefix",
        help="Output prefix for deduplicated fastq file(s).",
    ),
    statistics: str = typer.Option(
        "filtering_stats.json",
        "--statistics",
        help="Output path for stats file.",
    ),
    sample_name: str | None = typer.Option(
        None,
        "--sample-name",
        help="Name of sample e.g. DOX_treated_1.",
    ),
    read_type: str = typer.Option(
        "flashed",
        "--read-type",
        help="Type of read.",
    ),
    fragments: bool = typer.Option(
        True,
        "--fragments/--no-fragments",
        help="Determines if read fragment aggregations are produced.",
    ),
):
    """Remove unwanted aligned slices and identify reporters."""

    if method not in {"capture", "tri", "tiled"}:
        raise typer.BadParameter("method must be one of: capture, tri, tiled")

    if read_type.lower() not in {"flashed", "pe"}:
        raise typer.BadParameter("read-type must be one of: flashed, pe")

    from capcruncher.api.alignments_filter import filter as filter_slices

    filter_slices(
        method=method,
        bam=bam,
        annotations=annotations,
        custom_filtering=custom_filtering,
        output_prefix=output_prefix,
        statistics=statistics,
        sample_name=sample_name,
        read_type=read_type,
        fragments=fragments,
    )


cli = typer.main.get_command(alignments_app)
