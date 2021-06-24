import click


@click.group()
def cli():
    """
    Contains methods for genome digestion.
    """


@cli.command()
@click.argument("input_fasta")
@click.option(
    "-r", "--recognition_site", help="Recognition enzyme or sequence", required=True
)
@click.option(
    "-l", "--logfile", help="Path for digestion log file", default="genome_digest.log"
)
@click.option(
    "-o", "--output_file", help="Output file path", default="genome_digested.bed"
)
@click.option(
    "--remove_cutsite",
    help="Exclude the recognition sequence from the output",
    default=True,
)
@click.option(
    "--sort",
    help="Sorts the output bed file by chromosome and start coord.",
    default=False,
    is_flag=True,
)
def digest(*args, **kwargs):
    """
    Performs in silico digestion of a genome in fasta format.

    Digests the supplied genome fasta file and generates a bed file containing the
    locations of all restriction fragments produced by the supplied restriction enzyme.

    A log file recording the number of restriction fragments for the suplied genome is also
    generated.
    """
    from capcruncher.cli.genome_digest import digest

    digest(*args, **kwargs)