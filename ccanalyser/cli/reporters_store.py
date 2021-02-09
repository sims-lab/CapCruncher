import os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import click

from ccanalyser.cli import cli
from ccanalyser.tools.storage import create_cooler_cc, CoolerBinner


@cli.command()
@click.argument("counts")
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
    "-n", "--name", "capture_name", help="Name of capture oligo to store", required=True
)
@click.option(
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
def reporters_store(
    counts: str,
    restriction_fragment_map: str,
    capture_oligos: str,
    capture_name: str,
    binsizes: list = None,
    normalise=True,
    output="",
    overlap_fraction=0.5,
    genome: str = "",
    n_cores=2,
):
    """Stores interaction counts. Can also be used to bin interactions into constant genomic intervals"""

    # Load counts
    df_counts = pd.read_csv(counts, sep="\t")

    # Load restriction fragments
    df_restriction_fragment_map = pd.read_csv(
        restriction_fragment_map,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
    )

    # Create cooler file at restriction fragment resolution
    cooler_rf = create_cooler_cc(
        output,
        bins=df_restriction_fragment_map,
        pixels=df_counts,
        capture_name=capture_name,
        capture_oligos=capture_oligos,
    )

    # Bin cooler if required
    if binsizes:
        coolers_binned = [
            CoolerBinner(cooler_rf, binsize=b, n_cores=n_cores // len(binsizes))
            for b in binsizes
        ]
        Parallel(n_jobs=n_cores)(
            delayed(
                lambda binner: binner.to_cooler(
                    f"{output}_binned.hdf5", normalise=normalise
                )
            )
            for binner in coolers_binned
        )
