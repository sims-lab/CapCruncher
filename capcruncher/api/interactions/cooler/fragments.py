from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from capcruncher.api.interactions.cooler.create import create_cooler_cc


def fragments(
    counts: Path | str,
    fragment_map: Path | str,
    output: Path | str,
    viewpoint_path: Path | str,
    viewpoint_name: str = "",
    genome: str = "",
    suffix: str = "",
) -> None:
    """
    Store restriction-fragment interaction combinations in a cooler group.

    Parses reporter interaction counts and creates CapCruncher cooler output at
    restriction fragment resolution.
    """
    counts = os.fspath(counts)

    df_restriction_fragment_map = pd.read_csv(
        fragment_map,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
        dtype={"chrom": str},
    )

    if counts.endswith(".hdf5"):
        with pd.HDFStore(counts) as store:
            if not viewpoint_name:
                viewpoints = {key.split("/")[1] for key in store.keys()}
            else:
                viewpoints = {viewpoint_name}

            for viewpoint in viewpoints:
                create_cooler_cc(
                    output,
                    bins=df_restriction_fragment_map,
                    pixels=store[viewpoint],
                    viewpoint_name=viewpoint,
                    viewpoint_path=viewpoint_path,
                    assembly=genome,
                    suffix=suffix,
                )
        return

    create_cooler_cc(
        output,
        bins=df_restriction_fragment_map,
        pixels=pd.read_csv(counts, sep="\t"),
        viewpoint_name=viewpoint_name,
        viewpoint_path=viewpoint_path,
        assembly=genome,
        suffix=suffix,
    )
