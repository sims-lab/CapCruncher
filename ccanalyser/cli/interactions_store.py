import os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import click
import re

import cooler
import h5py
from ccanalyser.cli import cli
from ccanalyser.tools.storage import GenomicBinner, create_cooler_cc, CoolerBinner


@cli.group()
def interactions_store():
    """Stores interaction counts. Can also be used to bin interactions into constant genomic intervals"""


@interactions_store.command()
@click.argument("counts", required=True)
@click.option(
    "-f",
    "--fragment_map",
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
    "-n",
    "--capture_name",
    "capture_name",
    help="Name of capture oligo to store",
    required=True,
)
@click.option(
    "-g",
    "--genome",
    help="Name of genome",
)
@click.option(
    "--suffix",
    help="Suffix to append after the capture name for the output file",
)
@click.option(
    "-o",
    "--output",
    help="Name of output file. (Cooler formatted hdf5 file)",
    default="out.hdf5",
)
def fragments(
    counts, fragment_map, output, capture_name, capture_oligos, genome="", suffix=""
):

    # Load restriction fragments
    df_restriction_fragment_map = pd.read_csv(
        fragment_map,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
    )

    # Load counts
    df_counts = pd.read_csv(counts, sep="\t")

    # Create cooler file at restriction fragment resolution
    cooler_fn = create_cooler_cc(
        output,
        bins=df_restriction_fragment_map,
        pixels=df_counts,
        capture_name=capture_name,
        capture_oligos=capture_oligos,
        assembly=genome,
        suffix=suffix,
    )


@interactions_store.command()
@click.argument("cooler_fn", required=True)
@click.option(
    "-b",
    "--binsizes",
    help="Binsizes to use for windowing",
    default=(5000,),
    multiple=True,
    type=click.INT,
)
@click.option(
    "--normalise",
    is_flag=True,
    help="Enables normalisation of interaction counts during windowing",
)
@click.option(
    "--overlap_fraction",
    help="Minimum overlap between genomic bins and restriction fragments for overlap",
    default=0.5,
)
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
@click.option(
    "-o",
    "--output",
    help="Name of output file. (Cooler formatted hdf5 file)",
    default="out.hdf5",
)
def bins(
    cooler_fn,
    output,
    binsizes=None,
    normalise=False,
    n_cores=1,
    scale_factor=1e6,
    overlap_fraction=1e-9,
):

    
    for binsize in binsizes:
        cb = CoolerBinner(cooler_fn, binsize=binsize, n_cores=n_cores)
        cb.to_cooler(output, normalise=normalise, scale_factor=scale_factor)



@interactions_store.command()
@click.argument("coolers", required=True, nargs=-1)
@click.option("-o", "--output", help="Output file name")
def merge(coolers, output):


    with h5py.File(output, "w") as dest:

        for clr in coolers:
            re_fn = re.match("(.*)\.(.*)\.(.*)?\.hdf5", clr)
            sample = re_fn.group(1)
            capture = re_fn.group(2)
            resolution = re_fn.group(3)

            with h5py.File(clr, "r") as src:

                if resolution == "fragments":  # i.e. Is a fragment cooler
                    dest_grp_name = capture
                else:
                    dest_grp_name = f"{capture}/resolutions/{resolution}"

                if not dest.get(dest_grp_name):
                    dest.copy(src.parent, dest_grp_name)
                else:
                    for key in src.keys():
                        dest.copy(src[key], f'{dest_grp_name}/{key}')
                
                attributes = {k: v for k, v in src.parent.attrs.items()}
                dest[dest_grp_name].attrs.update(attributes)


    link_bins(output)


def link_bins(clr):
    """ Reduces cooler storage space by linking "bins" table within each hdf5 file as these are identical for all resolutions"""

    with h5py.File(clr, "a") as f:

        # Get all captures stored
        captures = list(f.keys())

        # Get all resolutions stored
        resolutions_group = f[captures[0]].get("resolutions")
        resolutions = list(resolutions_group.keys()) if resolutions_group else None

        for capture in captures[1:]:

            # Delete currenly stored bins group and replace with link to first capture "bins" group
            del f[capture]["bins"]
            f[capture]["bins"] = f[captures[0]]["bins"]

            if resolutions:
                for resolution in resolutions:
                    # Repeat for resolutions i.e. binned coolers
                    del f[capture]["resolutions"][resolution]["bins"]
                    f[capture]["resolutions"][resolution]["bins"] = f[captures[0]][
                        "resolutions"
                    ][resolution]["bins"]
