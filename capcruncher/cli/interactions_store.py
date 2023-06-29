from loguru import logger
import os
import tempfile
from typing import Tuple
import pandas as pd
import ray
from capcruncher.api.storage import (
    CoolerBinner,
    create_cooler_cc,
    merge_coolers,
)
import cooler
import warnings

# warnings.filterwarnings("ignore", category=DeprecationWarning)
# warnings.filterwarnings("ignore", category=FutureWarning)


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


@ray.remote(num_cpus=1)
def _bin_cooler(clr_in: os.PathLike, clr_out: os.PathLike, binsize: int, **kwargs):

    clr_binner = CoolerBinner(
        cooler_group=clr_in,
        binsize=binsize,
        **kwargs,
    )
    clr_binner.to_cooler(clr_out)
    return clr_out


def bins(
    cooler_path: os.PathLike,
    output: os.PathLike,
    binsizes: Tuple = None,
    normalise: bool = True,
    scale_factor: int = 1e6,
    overlap_fraction: float = 1e-9,
    conversion_tables: os.PathLike = None,
    n_cores: int = 1,
    **kwargs,
):
    """
    Convert a cooler group containing restriction fragments to constant genomic windows

    Parses a cooler group and aggregates restriction fragment interaction counts into
    genomic bins of a specified size. If the normalise option is selected,
    columns containing normalised counts are added to the pixels table of the output

    \f
    Args:
        cooler_path (os.PathLike): Path to cooler file.
        output (os.PathLike): Path to output cooler file.
        binsizes (Tuple, optional): Binsizes to bin cooler file to. Defaults to None.
        normalise (bool, optional): Whether to normalise counts. Defaults to True.
        scale_factor (int, optional): Scale factor for normalisation. Defaults to 1e6.
        overlap_fraction (float, optional): Minimum overlap fraction for binning. Defaults to 1e-9.
        conversion_tables (os.PathLike, optional): Path to conversion tables. Defaults to None.
        n_cores (int, optional): Number of cores to use. Defaults to 1.

    """
    clr_groups = cooler.api.list_coolers(cooler_path)

    assert clr_groups, "No cooler groups found in file"
    assert binsizes, "No binsizes provided"

    ray.init(num_cpus=n_cores, ignore_reinit_error=True)
    clr_tempfiles = []

    for binsize in binsizes:
        for clr_group in clr_groups:

            logger.info(f"Processing {clr_group}")
            clr_in = cooler.Cooler(f"{cooler_path}::{clr_group}")
            clr_out = tempfile.NamedTemporaryFile().name

            # TODO: Integrate these ino the CLI
            default_kwargs = dict(
                method="midpoint",
                minimum_overlap=0.51,
                n_cis_interaction_correction=True,
                n_rf_per_bin_correction=True,
                scale_factor=1_000_000,
            )

            clr_tempfiles.append(
                _bin_cooler.remote(clr_in, clr_out, binsize, **default_kwargs)
            )

    # Final cooler output
    clr_tempfiles = ray.get(clr_tempfiles)
    merge_coolers(clr_tempfiles, output)
