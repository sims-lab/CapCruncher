import typer

from capcruncher.cli.common import HELP_SETTINGS

genome_app = typer.Typer(
    help="Contains methods for genome digestion.",
    context_settings=HELP_SETTINGS,
    no_args_is_help=True,
)


@genome_app.callback()
def genome():
    """Contains methods for genome digestion."""


@genome_app.command()
def digest(
    input_fasta: str = typer.Argument(...),
    recognition_site: str = typer.Option(
        ...,
        "-r",
        "--recognition-site",
        "--recognition_site",
        help="Recognition enzyme or sequence.",
    ),
    logfile: str = typer.Option(
        "genome_digest.log",
        "-l",
        "--logfile",
        help="Path for digestion log file.",
    ),
    output_file: str = typer.Option(
        "genome_digested.bed",
        "-o",
        "--output-file",
        "--output_file",
        help="Output file path.",
    ),
    remove_cutsite: bool = typer.Option(
        True,
        "--remove-cutsite/--keep-cutsite",
        "--remove_cutsite/--keep_cutsite",
        help="Exclude the recognition sequence from the output.",
    ),
    sort: bool = typer.Option(
        False,
        "--sort",
        help="Sorts the output bed file by chromosome and start coord.",
    ),
):
    """
    Performs in silico digestion of a genome in fasta format.

    Digests the supplied genome fasta file and generates a bed file containing the
    locations of all restriction fragments produced by the supplied restriction enzyme.

    A log file recording the number of restriction fragments for the suplied genome is also
    generated.
    """
    from capcruncher.api.genome import digest_genome

    digest_genome(
        input_fasta=input_fasta,
        recognition_site=recognition_site,
        output_file=output_file,
        logfile=logfile,
        remove_cutsite=remove_cutsite,
        sort=sort,
    )


cli = typer.main.get_command(genome_app)
