import click

@click.group()
def cli():
    """Contains methods for fastq splitting, deduplicating and digestion."""


import ccanalyser.cli.fastq_split
import ccanalyser.cli.fastq_deduplicate
import ccanalyser.cli.fastq_digest


