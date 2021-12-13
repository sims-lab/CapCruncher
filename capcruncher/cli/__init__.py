import click
import os
from functools import cached_property
from importlib import import_module, metadata
import subprocess
import warnings
import logging


# create logger
logger = logging.getLogger("capcruncher")
logger.setLevel(logging.INFO)
logging.basicConfig(format='%(levelname)s:%(asctime)s %(module)-20s %(message)s', level=logging.INFO)

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


@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option("-h", "--help", is_flag=True)
@click.version_option(metadata.version(distribution_name="capcruncher"))
@click.argument("mode", type=click.Choice(["make", "run", "plot", "show", "clone", "touch"]))
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def pipeline(mode, pipeline_options, help=False, version=False):

    """Runs the data processing pipeline"""

    fn = os.path.abspath(__file__)
    dir_cli = os.path.dirname(fn)
    dir_package = os.path.dirname(dir_cli)

    cmd = [
        "python",
        f"{dir_package}/pipeline/pipeline.py",
        mode.replace("run", "make"),
    ]

    if help:
        cmd.append("--help")

    if pipeline_options:
        cmd.extend(pipeline_options)

    # Implicitly deal with the missing --local option
    if (
        not os.path.exists(os.environ.get("DRMAA_LIBRARY_PATH", ""))
        and not "--local" in pipeline_options
    ):
        warnings.showwarning(
            "DRMAA_LIBRARY_PATH is incorrect. Implicitly using --local with 4 cores",
            category=UserWarning,
            filename='CapCruncher CLI',
            lineno=113)
        cmd.append("--local")
        cmd.append("-p 4")

    completed = subprocess.run(cmd)

    if not completed.returncode == 0:
        raise RuntimeError("CapCruncher pipeline failed. Check pipeline.log for details")



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


@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_reporters:cli")
def reporters():
    """Reporter counting, storing, comparison and pileups"""

@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_plot:cli")
def plot():
    """
    Generates plots for the outputs produced by CapCruncher
    """

@cli.group(cls=LazyGroup, import_name="capcruncher.cli.cli_utilities:cli")
def utilities():
    """Contains miscellaneous functions"""

