import click

@click.group()
def cli():
    """Contains methods for reporter annotating, identifying and deduplication.
    """

import ccanalyser.cli.alignments_annotate
import ccanalyser.cli.alignments_deduplicate
import ccanalyser.cli.alignments_filter