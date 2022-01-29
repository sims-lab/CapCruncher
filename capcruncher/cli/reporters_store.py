import logging
import os
import pickle
from typing import Tuple

import h5py
import pandas as pd
from capcruncher.tools.storage import (
    CoolerBinner,
    create_cooler_cc,
    link_common_cooler_tables,
)


def get_viewpoints(store: h5py.File):

    keys = set()
    for k in store.keys():
        if "/" in k:
            keys.add(k.split("/")[1])
        else:
            keys.add(k)
    return keys


def get_dataset_keys(f):
    keys = []
    f.visit(lambda key: keys.append(key) if isinstance(f[key], h5py.Dataset) else None)
    return keys


def fragments(
    counts: os.PathLike,
    fragment_map: os.PathLike,
    output: os.PathLike,
    viewpoint_path: os.PathLike,
    viewpoint_name: str = "",
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

        with pd.HDFStore(counts) as store:

            if not viewpoint_name:
                viewpoints = {k.split("/")[1] for k in store.keys()}
            else:
                viewpoints = {
                    viewpoint_name,
                }

            for viewpoint in viewpoints:
                df_counts = store[viewpoint]

                create_cooler_cc(
                    output,
                    bins=df_restriction_fragment_map,
                    pixels=df_counts,
                    viewpoint_name=viewpoint,
                    viewpoint_path=viewpoint_path,
                    assembly=genome,
                    suffix=suffix,
                )

    else:
        df_counts = pd.read_csv(counts, sep="\t")
        # Create cooler file at restriction fragment resolution
        create_cooler_cc(
            output,
            bins=df_restriction_fragment_map,
            pixels=df_counts,
            viewpoint_name=viewpoint_name,
            viewpoint_path=viewpoint_path,
            assembly=genome,
            suffix=suffix,
        )


def bins(
    cooler_path: os.PathLike,
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

    # Get a conversion table if one exists
    if conversion_tables:
        with open(conversion_tables, "rb") as r:
            genomic_binner_objs = pickle.load(r)
    else:
        genomic_binner_objs = None

    # Get Viewpoints
    with h5py.File(cooler_path) as store:
        viewpoints = get_viewpoints(store)

    # Perform binning
    for viewpoint in viewpoints:

        cooler_group = f"{cooler_path}::/{viewpoint}"

        for binsize in binsizes:
            if genomic_binner_objs and (binsize in genomic_binner_objs):
                cb = CoolerBinner(
                    cooler_group,
                    binsize=binsize,
                    n_cores=n_cores,
                    binner=genomic_binner_objs[binsize],
                )
            else:
                cb = CoolerBinner(cooler_group, binsize=binsize, n_cores=n_cores)

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
    from collections import defaultdict

    import cooler

    logging.info("Merging cooler files")

    coolers_to_merge = defaultdict(list)

    # Remove output file as need to append to it.
    if os.path.exists(output):
        os.unlink(output)

    # Extract a list of coolers to merge, grouped by viewpoint name
    for clr in coolers:
        with h5py.File(clr, mode="r") as src:
            viewpoints = list(src.keys())

            for viewpoint in viewpoints:
                if not "resolutions" in list(src[viewpoint].keys()):
                    coolers_to_merge[viewpoint].append(f"{clr}::/{viewpoint}")
                else:
                    for resolution in src[viewpoint]["resolutions"].keys():
                        coolers_to_merge[f"{viewpoint}::{resolution}"].append(
                            f"{clr}::/{viewpoint}/resolutions/{resolution}"
                        )

    # Initial pass to perform copying for all coolers without a matching group
    need_merging = list()
    with h5py.File(output, mode="w") as dest:
        for ii, (viewpoint, cooler_uris) in enumerate(coolers_to_merge.items()):

            if len(cooler_uris) < 2:  # Only merge if two or more, else just copy
                (file_path, group_path) = cooler_uris[0].split("::")

                with h5py.File(file_path, mode="r") as src:
                    src.copy(src[group_path], dest, group_path)

            else:
                need_merging.append(viewpoint)

    # Actually merge the coolers left over that do have duplicates
    for viewpoint in need_merging:
        cooler_uris = coolers_to_merge[viewpoint]
        cooler.merge_coolers(
            f"{output}::/{viewpoint.replace('::', '/resolutions/')}", cooler_uris
        )

    # Reduce space by linking common tables (bins, chroms)
    link_common_cooler_tables(output)
