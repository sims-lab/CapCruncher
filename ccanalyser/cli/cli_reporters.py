import click 

@click.group()
def cli():
    """Contains methods for interaction counting, storing, bedgraph generation, comparisons and plotting.
    """

import ccanalyser.cli.reporters_pileup
import ccanalyser.cli.reporters_count
import ccanalyser.cli.reporters_differential
import ccanalyser.cli.reporters_heatmap
import ccanalyser.cli.reporters_store