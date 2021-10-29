import os
import re
from typing import Tuple
import warnings
warnings.simplefilter('ignore', category=RuntimeWarning)
import pickle

import click
import h5py
import pandas as pd
from capcruncher.cli.cli_reporters import cli
from capcruncher.tools.storage import CoolerBinner, create_cooler_cc, link_bins


def fragments(
    counts: os.PathLike,
    fragment_map: os.PathLike,
    output: os.PathLike,
    viewpoint_name: str,
    viewpoint_path: os.PathLike,
    genome: str = "",
    suffix: str = "",
):
    """
    Stores restriction fragment interaction combinations at the restriction fragment level.

    Parses reporter restriction fragment interaction counts produced by 
    "capcruncher reporters count" and gerates a cooler formatted group in an HDF5 File. 
    See `https://cooler.readthedocs.io/en/latest/` for further details.


    \f
    Args:
     counts (os.PathLike): Path to restriction fragment interactions counts .tsv file.
     fragment_map (os.PathLike): Path to restriction fragment .bed file, generated with genome-digest command.
     output (os.PathLike): Output file path for cooler hdf5 file.
     viewpoint_name (str): Name of viewpoint.
     viewpoint_path (os.PathLike): Path to viewpoints bed file.
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
    if counts.endswith(".hdf5"):
        df_counts = pd.read_hdf(counts, key=viewpoint_name)
    else:
        df_counts = pd.read_csv(counts, sep="\t")

    # Create cooler file at restriction fragment resolution
    cooler_fn = create_cooler_cc(
        output,
        bins=df_restriction_fragment_map,
        pixels=df_counts,
        viewpoint_name=viewpoint_name,
        viewpoint_path=viewpoint_path,
        assembly=genome,
        suffix=suffix,
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
     `:class:capcruncher.tools.storage.GenomicBinner` objects, one per binsize) can be supplied.


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
    else:
        genomic_binner_objs = None

    for binsize in binsizes:

        if genomic_binner_objs and (binsize in genomic_binner_objs):
            cb = CoolerBinner(cooler_fn, binsize=binsize, n_cores=n_cores, binner=genomic_binner_objs[binsize])
        else:
            cb = CoolerBinner(cooler_fn, binsize=binsize, n_cores=n_cores)

        if normalise:
            cb.normalise(scale_factor=scale_factor)

        cb.to_cooler(output)


def merge(coolers: Tuple, output: os.PathLike):
    """
    Merges capcruncher cooler files together.
    
    Produces a unified cooler with both restriction fragment and genomic bins whilst
    reducing the storage space required by hard linking the "bins" tables to prevent duplication.

    Args:
     coolers (Tuple): Cooler files produced by either the fragments or bins subcommands.
     output (os.PathLike): Path from merged cooler file.
    """    

    with h5py.File(output, "w") as dest:

        for clr in coolers:
            re_fn = re.match(r"(.*/)?(.*)\.(.*)\.(.*)?\.hdf5", clr)
            assert re_fn, f'{clr} file name not in correct format! Use format PATH_TO_HDF5/SAMPLE.VIEWPOINT.FRAGMENT|BINSIZE.hdf5'
            sample = re_fn.group(2)
            capture = re_fn.group(3)
            resolution = re_fn.group(4)

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


