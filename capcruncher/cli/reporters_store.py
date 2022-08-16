import logging
import os
import pickle
import tempfile
from typing import Iterable, Tuple

import h5py
import pandas as pd
from capcruncher.tools.storage import (
    CoolerBinner,
    GenomicBinner,
    create_cooler_cc,
    link_common_cooler_tables,
)
import cooler
import ujson

def get_merged_metadata(coolers: Iterable[os.PathLike]):
    """
    Merges metadata from multiple coolers.
    """
    # Get metadata from all coolers and copy to the merged file
    metadata = {}
    for cooler_uri in coolers:
        filepath, group = cooler_uri.split("::")

        with h5py.File(filepath, mode="r") as src:
            metadata_src = ujson.decode(src[group].attrs["metadata"])

            for metadata_key, metadata_value in metadata_src.items():

                if isinstance(metadata_value, str):
                    metadata[metadata_key] = metadata_value

                elif isinstance(metadata_value, Iterable):
                    if metadata_key not in metadata:
                        metadata[metadata_key] = []
                        metadata[metadata_key].extend(metadata_value)
                    else:
                        metadata[metadata_key].extend(
                            [
                                v
                                for v in metadata_value
                                if v not in metadata[metadata_key]
                            ]
                        )

                elif isinstance(metadata_value, (int, float)):
                    if metadata_key not in metadata:
                        metadata[metadata_key] = metadata_value
                    else:
                        metadata[metadata_key] += metadata_value

    return metadata


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
    normalise: bool = True,
    scale_factor: int = 1e6,
    overlap_fraction: float = 1e-9,
    conversion_tables: os.PathLike = None,
    **kwargs,
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
    clr_groups = cooler.api.list_coolers(cooler_path)

    if conversion_tables:
        logging.info("Loading conversion tables")
        with open(conversion_tables, "rb") as f:
            try:
                conversion_tables = pickle.load(f)
            except:
                raise ValueError("Conversion tables are not in a valid format.")
    else:
        logging.info("Generating conversion tables")
        clr_example = cooler.Cooler(f"{cooler_path}::{clr_groups[0]}")
        conversion_tables = {
            binsize: GenomicBinner(
                chromsizes=clr_example.chromsizes,
                fragments=clr_example.bins()[:],
                binsize=binsize,
            )
            for binsize in binsizes
        }

    clr_tempfiles = []
    for binsize in binsizes:
        for clr_group in clr_groups:

            binning_output = tempfile.NamedTemporaryFile().name

            logging.info(f"Processing {clr_group}")
            clr = cooler.Cooler(f"{cooler_path}::{clr_group}")
            clr_binner = CoolerBinner(
                cooler_group=clr, binner=conversion_tables[binsize]
            )

            if normalise:
                clr_binner.normalise(scale_factor=scale_factor)

            clr_binner.to_cooler(binning_output)
            clr_tempfiles.append(binning_output)

    # Final cooler output
    merge(clr_tempfiles, output)


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
        tmp = tempfile.NamedTemporaryFile().name
        cooler_uris = coolers_to_merge[viewpoint]
        cooler.merge_coolers(
            f"{tmp}::/{viewpoint.replace('::', '/resolutions/')}",
            cooler_uris,
            mergebuf=int(1e6),
        )

        with h5py.File(tmp, mode="r") as src:
            with h5py.File(output, mode="a") as dest:
                dest.copy(
                    src[viewpoint.replace("::", "/resolutions/")], dest, viewpoint
                )

        metadata = get_merged_metadata(cooler_uris)

        with h5py.File(output, mode="a") as dest:
            dest[viewpoint.replace("::", "/resolutions/")].attrs[
                "metadata"
            ] = ujson.encode(metadata)

    # Reduce space by linking common tables (bins, chroms)
    link_common_cooler_tables(output)
