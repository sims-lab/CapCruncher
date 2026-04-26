import click
from importlib import metadata

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def get_capcruncher_version() -> str:
    try:
        return metadata.version(distribution_name="capcruncher")
    except metadata.PackageNotFoundError:
        return "0+unknown"


class UnsortedGroup(click.Group):
    def list_commands(self, ctx):
        return list(self.commands)


@click.group(cls=UnsortedGroup)
@click.version_option(get_capcruncher_version())
def cli():
    """
    An end to end solution for processing: Capture-C, Tri-C and Tiled-C data.
    """


from capcruncher.cli.cli_alignments import cli as alignments_cli  # noqa: E402
from capcruncher.cli.cli_fastq import cli as fastq_cli  # noqa: E402
from capcruncher.cli.cli_genome import cli as genome_cli  # noqa: E402
from capcruncher.cli.cli_interactions import cli as interactions_cli  # noqa: E402
from capcruncher.cli.cli_plot import cli as plot_cli  # noqa: E402
from capcruncher.cli.cli_utilities import cli as utilities_cli  # noqa: E402

cli.add_command(fastq_cli, "fastq")
cli.add_command(genome_cli, "genome")
cli.add_command(alignments_cli, "alignments")
cli.add_command(interactions_cli, "interactions")
cli.add_command(plot_cli, "plot")
cli.add_command(utilities_cli, "utilities")


# Finally, import the pipeline command from the pipeline module
import capcruncher.cli.cli_pipeline  # noqa: E402,F401


__all__ = [
    "CONTEXT_SETTINGS",
    "UnsortedGroup",
    "cli",
    "get_capcruncher_version",
]
