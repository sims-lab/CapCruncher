import click

@click.group()
def cli():
    """Contains methods for reporter annotating, identifying and deduplication.
    """

import ccanalyser.cli.reporters_annotate
import ccanalyser.cli.reporters_deduplicate
import ccanalyser.cli.reporters_identify