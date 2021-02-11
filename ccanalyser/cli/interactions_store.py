import os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import click

from ccanalyser.cli import cli
from ccanalyser.tools.storage import create_cooler_cc, CoolerBinner


@cli.command()
@click.argument("counts_fnames", nargs=-1, required=True)
@click.option(
    "-r",
    "--restriction_fragment_map",
    help="Path to digested genome bed file",
    required=True,
)
@click.option(
    "-c",
    "--capture_oligos",
    "capture_oligos",
    help="Path to capture oligos file",
    required=True,
)
@click.option(
    "-n", "--names", "capture_names", help="Name of capture oligo to store", required=True, multiple=True,
)
@click.option(
    '-b',
    "--binsizes",
    help="Binsizes to use for windowing",
    default=0,
    multiple=True,
    type=click.INT,
)
@click.option(
    "-o",
    "--output",
    help="Name of output file. (Cooler formatted hdf5 file)",
    default="out.hdf5",
)
@click.option('--normalise', is_flag=True, help='Enables normalisation of interaction counts during windowing')
@click.option(
    "--overlap_fraction",
    help="Minimum overlap between genomic bins and restriction fragments for overlap",
    default=0.5,
)
@click.option("-g", "--genome", help="Genome used for alignment")
@click.option(
    "-p",
    "--n_cores",
    help="Number of cores used for binning",
    default=4,
    type=click.INT,
)
@click.option(
    "--scale_factor",
    help="Scaling factor used for normalisation",
    default=1e6,
    type=click.INT,
)
def interactions_store(
    counts_fnames: tuple,
    restriction_fragment_map: str,
    capture_oligos: str,
    capture_names: str,
    binsizes: list = None,
    normalise=False,
    output="",
    overlap_fraction=0.5,
    genome: str = "",
    n_cores=2,
    scale_factor=1e6,
):
    """Stores interaction counts. Can also be used to bin interactions into constant genomic intervals"""
   

    # Load restriction fragments
    df_restriction_fragment_map = pd.read_csv(
        restriction_fragment_map,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
    )


    for capture_name, counts in zip(capture_names, counts_fnames):


        # Load counts
        df_counts = pd.read_csv(counts, sep="\t")


        # Create cooler file at restriction fragment resolution
        cooler_rf = create_cooler_cc(
            output,
            bins=df_restriction_fragment_map,
            pixels=df_counts,
            capture_name=capture_name,
            capture_oligos=capture_oligos,
            assembly=genome
        )


        if binsizes:

            for b in binsizes:
                binner = CoolerBinner(cooler_rf, binsize=b, n_cores=n_cores, normalise=normalise, scale_factor=scale_factor)
                binner.to_cooler(output, normalise=normalise)
        


