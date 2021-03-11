import click
import os
from functools import cached_property
from pkg_resources import iter_entry_points
from importlib import import_module
import subprocess


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
        module, name = self._import_name.split(':', 1)
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
def cli():
    """

    An end to end solution for processing:
    
    Capture-C, Tri-C and Tiled-C data.

    Type -h or --help after any subcommand for more information.
    """


@cli.group(cls=LazyGroup, import_name='ccanalyser.cli._fastq:cli')
def fastq():
    """Fastq splitting, deduplicating and digestion.
    """

@cli.group(cls=LazyGroup, import_name='ccanalyser.cli._genome:cli')
def genome():
    """Genome digestion.
    """

@cli.group(cls=LazyGroup, import_name='ccanalyser.cli._reporters:cli')
def reporters():
    """Reporter annotation, identification and deduplication.
    """

@cli.group(cls=LazyGroup, import_name='ccanalyser.cli._interactions:cli')
def interactions():
    """Interaction counting, storing, comparison, plotting and bedgraph generation.
    """

@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option("-h", "--help", default=False, is_flag=True)
@click.argument("mode", type=click.Choice(["make", "show", "clone", "touch"]))
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def pipeline(mode, pipeline_options, help=False):
    '''Runs the data processing pipeline'''
    
    fn = os.path.abspath(__file__)
    dir_cli = os.path.dirname(fn)
    dir_package = os.path.dirname(dir_cli)
    

    cmd = ['python', 
           f'{dir_package}/pipeline/pipeline.py', 
           mode,
           ]
    
    if pipeline_options:
        cmd.extend(pipeline_options.split())
    
    subprocess.run(cmd)



# TODO: LAZY imports to speed up cli

# from . import (genome_digest,
#                fastq_split,
#                fastq_deduplicate, 
#                fastq_digest,
#                slices_annotate,
#                reporters_identify,
#                reporters_deduplicate,
#                interactions_count,
#                interactions_store,
#                interactions_bedgraph,
#                interactions_plot_dev,
#                interactions_differential,
#                )




    




