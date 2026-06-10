from __future__ import annotations

import json
import os
import tempfile
from collections.abc import Iterable
from pathlib import Path

import h5py
from loguru import logger


def link_common_cooler_tables(clr: Path | str) -> None:
    """Reduces cooler storage space by linking "bins" table.

     All of the cooler "bins" tables containing the genomic coordinates of each bin
     are identical for all cooler files of the same resoultion. As cooler.create_cooler
     generates a new bins table for each cooler, this leads to a high degree of duplication.

     This function hard links the bins tables for a given resolution to reduce the degree of duplication.

    Args:
     clr (os.PathLike): Path to cooler hdf5 produced by the merge command.
    """

    logger.info("Making links to common cooler tables to conserve disk space")

    with h5py.File(clr, "a") as f:
        # Get all viewpoints stored
        viewpoints = sorted(list(f.keys()))

        # Get all resolutions stored
        try:
            resolutions = [res for res in f[viewpoints[0]]["resolutions"]]
        except (KeyError, IndexError):
            resolutions = None

        for viewpoint in viewpoints[1:]:
            try:
                # Delete currenly stored bins group and replace with link to first viewpoint "bins" group
                del f[viewpoint]["bins"]
                f[viewpoint]["bins"] = f[viewpoints[0]]["bins"]

                # Delete chroms table and replace with link to the first "chroms" group
                del f[viewpoint]["chroms"]
                f[viewpoint]["chroms"] = f[viewpoints[0]]["chroms"]
            except KeyError:
                pass

            # Repeat for resolutions i.e. binned coolers
            if resolutions:
                for resolution in resolutions:
                    del f[viewpoint]["resolutions"][resolution]["bins"]
                    f[viewpoint]["resolutions"][resolution]["bins"] = f[viewpoints[0]][
                        "resolutions"
                    ][resolution]["bins"]

                    del f[viewpoint]["resolutions"][resolution]["chroms"]
                    f[viewpoint]["resolutions"][resolution]["chroms"] = f[
                        viewpoints[0]
                    ]["resolutions"][resolution]["chroms"]


def get_merged_cooler_metadata(coolers: Iterable[Path | str]) -> dict:
    """
    Merges metadata from multiple coolers.
    """
    # Get metadata from all coolers and copy to the merged file
    metadata = {}
    for cooler_uri in coolers:
        filepath, group = os.fspath(cooler_uri).split("::")

        with h5py.File(filepath, mode="r") as src:
            metadata_src = json.loads(src[group].attrs["metadata"])

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


def merge_coolers(
    coolers: tuple[Path | str, ...] | list[Path | str], output: Path | str
):
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

    logger.info("Merging cooler files")

    coolers_to_merge = defaultdict(list)

    # Remove output file as need to append to it.
    if os.path.exists(output):
        os.unlink(output)

    # Extract a list of coolers to merge, grouped by viewpoint name
    for clr in coolers:
        with h5py.File(clr, mode="r") as src:
            viewpoints = list(src.keys())

            for viewpoint in viewpoints:
                if "resolutions" not in list(src[viewpoint].keys()):
                    coolers_to_merge[viewpoint].append(f"{clr}::/{viewpoint}")
                else:
                    for resolution in src[viewpoint]["resolutions"].keys():
                        coolers_to_merge[f"{viewpoint}::{resolution}"].append(
                            f"{clr}::/{viewpoint}/resolutions/{resolution}"
                        )

    # Initial pass to perform copying for all coolers without a matching group
    need_merging = list()
    with h5py.File(output, mode="w") as dest:
        for _, (viewpoint, cooler_uris) in enumerate(coolers_to_merge.items()):
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

        metadata = get_merged_cooler_metadata(cooler_uris)

        with h5py.File(output, mode="a") as dest:
            dest[viewpoint.replace("::", "/resolutions/")].attrs["metadata"] = (
                json.dumps(metadata)
            )

    # Reduce space by linking common tables (bins, chroms)
    link_common_cooler_tables(output)
