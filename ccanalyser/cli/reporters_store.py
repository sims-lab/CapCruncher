import os
import re
from typing import Tuple
import warnings
warnings.simplefilter('ignore', category=RuntimeWarning)
import pickle

import click
import h5py
import pandas as pd
from ccanalyser.cli.cli_reporters import cli
from ccanalyser.tools.storage import CoolerBinner, create_cooler_cc, link_bins


@cli.group()
def store():
    """
    Store reporter counts.

    These commands store and manipulate reporter restriction fragment interaction 
    counts as cooler formated groups in HDF5 files.

    See subcommands for details. 

    """


@store.command()
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
    counts: os.PathLike,
    fragment_map: os.PathLike,
    output: os.PathLike,
    capture_name: str,
    capture_oligos: os.PathLike,
    genome: str = "",
    suffix: str = "",
):
    """
    Stores restriction fragment interaction combinations at the restriction fragment level.

    Parses reporter restriction fragment interaction counts produced by 
    "ccanalyser reporters count" and gerates a cooler formatted group in an HDF5 File. 
    See `https://cooler.readthedocs.io/en/latest/` for further details.


    \f
    Args:
     counts (os.PathLike): Path to restriction fragment interactions counts .tsv file.
     fragment_map (os.PathLike): Path to restriction fragment .bed file, generated with genome-digest command.
     output (os.PathLike): Output file path for cooler hdf5 file.
     capture_name (str): Name of capture probe.
     capture_oligos (os.PathLike): Path to capture oligos bed file.
     genome (str, optional): Name of genome used for alignment e.g. hg19. Defaults to "".
     suffix (str, optional): Suffix to append to filename. Defaults to "".
    """
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


@store.command()
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
    "--conversion_tables",
    help="Pickle file containing pre-computed fragment -> bin conversions.",
    default=None,
)
@click.option(
    "-o",
    "--output",
    help="Name of output file. (Cooler formatted hdf5 file)",
    default="out.hdf5",
)
def bins(
    cooler_fn: os.PathLike,
    output: os.PathLike,
    binsizes: Tuple = None,
    normalise: bool = False,
    n_cores: int = 1,
    scale_factor: int = 1e6,
    overlap_fraction: float = 1e-9,
    conversion_tables: os.PathLike = None,
):
    """
    Convert a cooler group containing restriction fragments to constant genomic windows

    Parses a cooler group and aggregates restriction fragment interaction counts into 
    genomic bins of a specified size. If the normalise option is selected,  
    columns containing normalised counts are added to the pixels table of the output 


    Notes:
     To avoid repeatedly calculating restriction fragment to bin conversions,
     bin conversion tables (a .pkl file containing a dictionary of 
     `:class:ccanalyser.tools.storage.GenomicBinner` objects, one per binsize) can be supplied.


    Args:
     cooler_fn (os.PathLike): Path to cooler file. Nested coolers can be specified by STORE_FN.hdf5::/PATH_TO_COOLER
     output (os.PathLike): Path for output binned cooler file.
     binsizes (Tuple, optional): Genomic window sizes to use for binning. Defaults to None.
     normalise (bool, optional): Normalise the number of interactions to total number of cis interactions (True). Defaults to False.
     n_cores (int, optional): Number of cores to use for binning. Performed in parallel by chromosome. Defaults to 1.
     scale_factor (int, optional): Scaling factor to use for normalising interactions. Defaults to 1e6.
     overlap_fraction (float, optional): Minimum fraction to use for defining overlapping bins. Defaults to 1e-9.
    """

    if conversion_tables:
        with open(conversion_tables, 'rb') as r:
            genomic_binner_objs = pickle.load(r)

    for binsize in binsizes:

        if binsize in genomic_binner_objs:
            cb = CoolerBinner(cooler_fn, binsize=binsize, n_cores=n_cores, binner=genomic_binner_objs[binsize])
        else:
            cb = CoolerBinner(cooler_fn, binsize=binsize, n_cores=n_cores)

        cb.to_cooler(output, normalise=normalise, scale_factor=scale_factor)


@store.command()
@click.argument("coolers", required=True, nargs=-1)
@click.option("-o", "--output", help="Output file name")
def merge(coolers: Tuple, output: os.PathLike):
    """
    Merges ccanalyser cooler files together.
    
    Produces a unified cooler with both restriction fragment and genomic bins whilst
    reducing the storage space required by hard linking the "bins" tables to prevent duplication.

    Args:
     coolers (Tuple): Cooler files produced by either the fragments or bins subcommands.
     output (os.PathLike): Path from merged cooler file.
    """    

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
                        dest.copy(src[key], f"{dest_grp_name}/{key}")

                attributes = {k: v for k, v in src.parent.attrs.items()}
                dest[dest_grp_name].attrs.update(attributes)

    link_bins(output)


