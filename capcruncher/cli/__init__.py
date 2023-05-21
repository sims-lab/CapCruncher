import click
from functools import cached_property
from importlib import import_module, metadata
from loguru import logger

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


class UnsortedGroup(click.Group):
    def list_commands(self, ctx):
        return list(self.commands)


class LazyGroup(click.Group):
    """
    A click Group that imports the actual implementation only when
    needed.  This allows for more resilient CLIs where the top-level
    command does not fail when a subcommand is broken enough to fail
    at import time.
    """

    def __init__(self, import_name, **kwargs):
        self._import_name = import_name
        super().__init__(**kwargs)

    @cached_property
    def _impl(self):
        module, name = self._import_name.split(":", 1)
        return getattr(import_module(module), name)

    def get_command(self, ctx, cmd_name):
        return self._impl.get_command(ctx, cmd_name)

    def list_commands(self, ctx):
        return self._impl.list_commands(ctx)

    def invoke(self, ctx):
        return self._impl.invoke(ctx)

    def get_usage(self, ctx):
        return self._impl.get_usage(ctx)

    def get_params(self, ctx):
        return self._impl.get_params(ctx)


@click.group(cls=UnsortedGroup)
@click.version_option(metadata.version(distribution_name="capcruncher"))
def cli():
    """
    An end to end solution for processing: Capture-C, Tri-C and Tiled-C data.
    """


@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_fastq:cli")
def fastq():
    """
    Fastq splitting, deduplication and digestion.
    """


@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_genome:cli")
def genome():
    """
    Genome wide methods digestion.
    """


@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_alignments:cli")
def alignments():
    """Alignment annotation, identification and deduplication."""


@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_interactions:cli")
def interactions():
    """Reporter counting, storing, comparison and pileups"""


@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_plot:cli")
def plot():
    """
    Generates plots for the outputs produced by CapCruncher
    """


@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_utilities:cli")
def utilities():
    """Contains miscellaneous functions"""


# Finally, import the pipeline command from the pipeline module
import capcruncher.cli.cli_pipeline


__all__ = [
    "alignments_annotate",
    "alignments_deduplicate",
    "alignments_filter",
    "fastq_deduplicate",
    "fastq_split",
    "fastq_digest",
    "fastq_split",
    "genome_digest",
    "plot",
    "reporters_compare",
    "reporters_count",
    "reporters_differential",
    "reporters_pileup",
    "reporters_store",
]
