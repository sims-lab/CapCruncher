from __future__ import annotations

import os
from pathlib import Path

import cooler
import pandas as pd
import pyranges1 as pr

from capcruncher.types import Assay
from capcruncher.api.interactions.cooler.viewpoints import Viewpoint

PIXEL_COLUMNS = ["bin1_id", "bin2_id", "count"]


def _normalise_pixels(pixels: pd.DataFrame) -> pd.DataFrame:
    pixels = pd.DataFrame(pixels).copy()

    if pixels.empty:
        raise ValueError("Cannot create a cooler group from an empty pixels table.")

    missing_columns = [
        column for column in PIXEL_COLUMNS if column not in pixels.columns
    ]
    if missing_columns:
        raise ValueError(
            "Pixels table is missing required column(s): "
            + ", ".join(missing_columns)
        )

    return pixels


def create_cooler_cc(
    output_prefix: Path | str,
    bins: pd.DataFrame,
    pixels: pd.DataFrame,
    viewpoint_name: str,
    viewpoint_path: Path | str,
    assay: Assay | str = Assay.CAPTURE,
    suffix: str | None = None,
    **cooler_kwargs,
) -> str:
    """
    Creates a cooler hdf5 file or cooler formatted group within a hdf5 file.

    Args:
     output_prefix (str): Output path for hdf5 file. If this already exists, will append a new group to the file.
     bins (pd.DataFrame): DataFrame containing the genomic coordinates of all bins in the pixels table.
     pixels (pd.DataFrame): DataFrame with columns: bin1_id, bin2_id, count.
     viewpoint_name (str): Name of viewpoint to store.
     viewpoint_path (os.PathLike): Path to viewpoints used for the analysis.
     suffix (str, optional): Suffix to append before the .hdf5 file extension. Defaults to None.

    Raises:
     ValueError: Viewpoint name must exactly match the a supplied viewpoint.

    Returns:
     os.PathLike: Path of cooler hdf5 file.
    """
    output_prefix = os.fspath(output_prefix)

    viewpoint = Viewpoint.from_bed(
        bed=viewpoint_path, viewpoint=viewpoint_name, assay=assay
    )

    bins = pd.DataFrame(bins).copy()
    pixels = _normalise_pixels(pixels)

    gr_bins = pr.PyRanges(
        bins.rename(
            columns={
                "chrom": "Chromosome",
                "start": "Start",
                "end": "End",
                "name": "Name",
            }
        )
    )

    # Get cis bins
    bins_cis = viewpoint.bins_cis(gr_bins)

    # Get cis pixels
    pixels_cis = pixels.loc[
        lambda df: (df["bin1_id"].isin(bins_cis)) | (df["bin2_id"].isin(bins_cis))
    ]

    # Metadata for cooler file.
    metadata = {
        "viewpoint_bins": viewpoint.bin_names(gr_bins),
        "viewpoint_name": viewpoint_name,
        "viewpoint_chrom": viewpoint.chromosomes,
        "viewpoint_coords": viewpoint.coords,
        "n_cis_interactions": int(pixels_cis["count"].sum()),
        "n_total_interactions": int(pixels["count"].sum()),
    }

    if os.path.exists(
        output_prefix
    ):  # Will append to a prexisting file if one is supplied
        append_to_file = True
        cooler_fn = f"{output_prefix}::/{viewpoint_name}"
    else:
        append_to_file = False
        cooler_fn = f"{output_prefix.replace('.hdf5', '')}{'.' + suffix if suffix else ''}.hdf5::/{viewpoint_name}"

    cooler.create_cooler(
        cooler_fn,
        bins=bins,
        pixels=pixels,
        metadata=metadata,
        mode="w" if not append_to_file else "a",
        **cooler_kwargs,
    )

    return cooler_fn
