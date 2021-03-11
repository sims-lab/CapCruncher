import click 

@click.group()
def cli():
    """Contains methods for interaction counting, storing, bedgraph generation, comparisons and plotting.
    """

import ccanalyser.cli.interactions_bedgraph
import ccanalyser.cli.interactions_count
import ccanalyser.cli.interactions_differential
import ccanalyser.cli.interactions_plot
import ccanalyser.cli.interactions_store