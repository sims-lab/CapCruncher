from importlib import metadata

import typer

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def get_capcruncher_version() -> str:
    try:
        return metadata.version(distribution_name="capcruncher")
    except metadata.PackageNotFoundError:
        return "0+unknown"


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(get_capcruncher_version())
        raise typer.Exit()


app = typer.Typer(
    help="An end to end solution for processing: Capture-C, Tri-C and Tiled-C data.",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
)


@app.callback()
def capcruncher(
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the version and exit.",
    ),
) -> None:
    """An end to end solution for processing: Capture-C, Tri-C and Tiled-C data."""


from capcruncher.cli.alignments import alignments_app  # noqa: E402
from capcruncher.cli.fastq import fastq_app  # noqa: E402
from capcruncher.cli.genome import genome_app  # noqa: E402
from capcruncher.cli.interactions import interactions_app  # noqa: E402
from capcruncher.cli.pipeline import pipeline_app, pipeline_config, pipeline_init  # noqa: E402
from capcruncher.cli.plot import plot_app  # noqa: E402
from capcruncher.cli.utilities import utilities_app  # noqa: E402

app.add_typer(fastq_app, name="fastq")
app.add_typer(genome_app, name="genome")
app.add_typer(alignments_app, name="alignments")
app.add_typer(interactions_app, name="interactions")
app.add_typer(pipeline_app, name="pipeline")
app.command(name="pipeline-init")(pipeline_init)
app.command(name="pipeline-config")(pipeline_config)
app.add_typer(plot_app, name="plot")
app.add_typer(utilities_app, name="utilities")

cli = typer.main.get_command(app)


__all__ = [
    "CONTEXT_SETTINGS",
    "app",
    "cli",
    "get_capcruncher_version",
]
