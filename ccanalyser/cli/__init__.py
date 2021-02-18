import click


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


class UnsortedGroup(click.Group):
    def list_commands(self, ctx):
        return list(self.commands)


@click.group(context_settings=CONTEXT_SETTINGS, cls=UnsortedGroup)
def cli():
    """
    Type -h or --help after any subcommand for more information.
    """
    pass

from . import (genome_digest,
               fastq_split,
               fastq_deduplicate, 
               fastq_digest,
               slices_annotate,
               reporters_identify,
               reporters_deduplicate,
               interactions_count,
               interactions_store,
               interactions_bedgraph,
               )

from ccanalyser.pipeline.pipeline import pipeline


