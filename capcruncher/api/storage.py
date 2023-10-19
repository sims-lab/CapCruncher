import os
import tempfile
import pandas as pd
import numpy as np
from pybedtools import BedTool
import cooler
import h5py
import functools
import itertools
from loguru import logger
import ujson
from typing import Iterable, Tuple, Union, List, Dict, Literal
import pyranges as pr
import re
from dataclasses import dataclass


class Viewpoint:
    def __init__(
        self, coordinates: pr.PyRanges, assay: Literal["capture", "tri", "tiled"]
    ) -> None:
        self.coordinates = coordinates
        self.assay = assay

    @classmethod
    def from_bed(
        cls, bed: str, viewpoint: str, assay: Literal["capture", "tri", "tiled"]
    ):
        """
        Creates a viewpoint object from a bed file.

        Args:
            bed (str): Path to bed file containing viewpoint coordinates.
            viewpoint (str): Name of viewpoint to extract from bed file.

        Raises:
            IndexError: Oligo name cannot be found within viewpoints.

        Returns:
            Viewpoint: Viewpoint object.
        """
        gr_viewpoints = pr.read_bed(bed)
        df_viewpoints = gr_viewpoints.as_df()

        df_viewpoints = df_viewpoints.loc[
            lambda df: df["Name"].str.contains(f"{viewpoint}$")
        ]

        if df_viewpoints.empty:
            raise IndexError(
                f"Oligo name cannot be found within viewpoints: {viewpoint}"
            )

        return Viewpoint(df_viewpoints.pipe(pr.PyRanges), assay=assay)

    def bins(self, bins: pr.PyRanges):
        """
        Returns the bins that overlap with the viewpoint.

        Args:
            bins (pr.PyRanges): PyRanges object containing all bins.

        Returns:
            pr.PyRanges: PyRanges object containing all bins that overlap with the viewpoint.
        """
        return bins.join(self.coordinates)

    def bin_names(self, bins: pr.PyRanges) -> List[int]:
        return self.bins(bins).df["Name"].astype(int).to_list()

    def bins_cis(self, bins: pr.PyRanges) -> List[int]:
        """
        Returns the bins that are on the same chromosome(s) as the viewpoint.

        Args:
            bins (pr.PyRanges): PyRanges object containing all bins.

        Returns:
            List[int]: List of bin names.
        """

        # Get the chromosomes of the viewpoint
        viewpoint_chromosomes = self.chromosomes

        # Get the bins that are on the same chromosome(s) as the viewpoint
        df_cis_bins = bins.df.loc[
            lambda df: df["Chromosome"].isin(viewpoint_chromosomes)
        ]

        # If capture or tri, remove viewpoint bins from cis bins
        if self.assay == "capture" or self.assay == "tri":
            df_cis_bins = df_cis_bins.loc[
                lambda df: ~df["Name"].isin(self.bin_names(bins))
            ]

        return df_cis_bins["Name"].to_list()

    @property
    def chromosomes(self) -> List[str]:
        return self.coordinates.df["Chromosome"].unique().tolist()

    @property
    def coords(self) -> List[str]:
        """
        Returns the genomic coordinates of the viewpoint.

        Returns:
            List[str]: List of genomic coordinates.
        """
        _coords = []
        for row in self.coordinates.df.itertuples():
            _coords.append(f"{row.Chromosome}:{row.Start}-{row.End}")

        return _coords


def create_cooler_cc(
    output_prefix: str,
    bins: pd.DataFrame,
    pixels: pd.DataFrame,
    viewpoint_name: str,
    viewpoint_path: os.PathLike,
    assay: Literal["capture", "tri", "tiled"] = "capture",
    suffix=None,
    **cooler_kwargs,
) -> os.PathLike:
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

    viewpoint = Viewpoint.from_bed(
        bed=viewpoint_path, viewpoint=viewpoint_name, assay=assay
    )

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


class CoolerBinner:
    def __init__(
        self,
        cooler_group: os.PathLike,
        binsize: int = None,
        method: Union[Literal["overlap"], Literal["midpoint"]] = "midpoint",
        minimum_overlap: float = 0.51,
        n_cis_interaction_correction: bool = True,
        n_rf_per_bin_correction: bool = True,
        scale_factor: int = 1_000_000,
        assay: Literal["capture", "tri", "tiled"] = "capture",
    ) -> None:
        self.cooler_group = cooler_group
        self.binsize = binsize
        self.method = method
        self.minimum_overlap = minimum_overlap

        if isinstance(cooler_group, str):
            self.cooler = cooler.Cooler(cooler_group)
        elif isinstance(cooler_group, cooler.Cooler):
            self.cooler = cooler_group
        else:
            raise ValueError(
                "cooler_group must be a path to a cooler file or a cooler object"
            )

        self.n_cis_interactions = self.cooler.info["metadata"]["n_cis_interactions"]
        self.n_cis_interaction_correction = n_cis_interaction_correction
        self.n_restriction_fragment_correction = n_rf_per_bin_correction
        self.scale_factor = scale_factor
        self.assay = assay

    @functools.cached_property
    def genomic_bins(self) -> pr.PyRanges:
        return (
            cooler.binnify(binsize=self.binsize, chromsizes=self.cooler.chromsizes)
            .sort_values(by=["chrom", "start", "end"])
            .assign(
                genomic_bin_id=lambda df: df.reset_index(drop=True)
                .index.to_series()
                .values
            )
            .rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
            .pipe(pr.PyRanges)
        )

    @functools.cached_property
    def fragment_bins(self):
        return (
            self.cooler.bins()[:]
            .rename(
                columns={
                    "chrom": "Chromosome",
                    "start": "Start",
                    "end": "End",
                    "name": "fragment_id",
                }
            )
            .pipe(pr.PyRanges)
        )

    @functools.cached_property
    def fragment_to_genomic_table(self) -> pr.PyRanges:
        """
        Translate genomic bins to fragment bins
        """

        fragment_bins = self.fragment_bins

        if self.method == "midpoint":
            fragment_bins = (
                fragment_bins.as_df()
                .assign(
                    Start=lambda df: df["Start"] + (df["End"] - df["Start"]) / 2,
                    End=lambda df: df["Start"] + 1,
                )
                .pipe(pr.PyRanges)
            )

        pr_fragment_to_bins = self.genomic_bins.join(
            fragment_bins, strandedness=0, how=None, report_overlap=True
        )

        if self.method == "overlap":
            pr_fragment_to_bins = pr_fragment_to_bins[
                pr_fragment_to_bins["Overlap"] >= self.minimum_overlap
            ]

        # Add number of fragments per bin
        pr_fragment_to_bins = pr_fragment_to_bins.assign(
            "n_fragments_per_bin",
            lambda df: df.groupby("genomic_bin_id")["fragment_id"].transform("nunique"),
        )

        return pr_fragment_to_bins

    @functools.cached_property
    def fragment_to_genomic_mapping(self) -> Dict[int, int]:
        """
        Translate genomic bins to fragment bins
        """
        fragment_to_bins_mapping = (
            self.fragment_to_genomic_table.as_df()
            .set_index("fragment_id")["genomic_bin_id"]
            .to_dict()
        )
        return fragment_to_bins_mapping

    @functools.cached_property
    def pixels(self) -> pd.DataFrame:
        """
        Translate fragment pixels to genomic pixels
        """

        fragment_to_bins_mapping = self.fragment_to_genomic_mapping

        pixels = self.cooler.pixels()[:].assign(
            genomic_bin1_id=lambda df: df["bin1_id"].map(fragment_to_bins_mapping),
            genomic_bin2_id=lambda df: df["bin2_id"].map(fragment_to_bins_mapping),
        )

        # Sum the counts of pixels that map to the same genomic bins
        pixels = (
            pixels.groupby(["genomic_bin1_id", "genomic_bin2_id"])
            .agg(
                count=("count", "sum"),
            )
            .reset_index()
        )

        # Normalize pixels if specified
        if self.n_restriction_fragment_correction:
            n_fragments_per_bin = (
                self.fragment_to_genomic_table.as_df()
                .set_index("genomic_bin_id")["n_fragments_per_bin"]
                .to_dict()
            )
            pixels = pixels.assign(
                n_fragments_per_bin1=lambda df: df["genomic_bin1_id"].map(
                    n_fragments_per_bin
                ),
                n_fragments_per_bin2=lambda df: df["genomic_bin2_id"].map(
                    n_fragments_per_bin
                ),
                n_fragments_per_bin_correction=lambda df: (
                    df["n_fragments_per_bin1"] + df["n_fragments_per_bin2"]
                ),
                count_n_rf_norm=lambda df: df["count"]
                / df["n_fragments_per_bin_correction"],
            )

        if self.n_cis_interaction_correction:
            pixels = pixels.assign(
                count_n_cis_norm=lambda df: (df["count"] / self.n_cis_interactions)
                * self.scale_factor,
            )

        if self.n_cis_interaction_correction and self.n_restriction_fragment_correction:
            pixels = pixels.assign(
                count_n_cis_rf_norm=lambda df: (
                    pixels["count_n_rf_norm"] / self.n_cis_interactions
                )
                * self.scale_factor
            )

        return pixels

    @functools.cached_property
    def viewpoint_bins(self) -> List[int]:
        """
        Return list of viewpoint bins
        """

        pr_viewpoint = pr.from_dict(
            dict(
                zip(
                    ["Chromosome", "Start", "End"],
                    [
                        [
                            x,
                        ]
                        for x in re.split(
                            ":|-", self.cooler.info["metadata"]["viewpoint_coords"][0]
                        )
                    ],
                )
            )
        )

        return pr_viewpoint.join(self.genomic_bins).df["genomic_bin_id"].to_list()

    def to_cooler(self, store: os.PathLike):
        metadata = {**self.cooler.info["metadata"]}
        metadata["viewpoint_bins"] = [int(x) for x in self.viewpoint_bins]
        metadata["n_interactions_total"] = int(self.cooler.pixels()[:]["count"].sum())
        cooler_fn = f"{store}::/{metadata['viewpoint_name']}/resolutions/{self.binsize}"

        pixels = (
            self.pixels.drop(
                columns=[
                    "bin1_id",
                    "bin2_id",
                    "n_fragments_per_bin1",
                    "n_fragments_per_bin2",
                    "n_fragments_per_bin_correction",
                ],
                errors="ignore",
            )
            .rename(
                columns={"genomic_bin1_id": "bin1_id", "genomic_bin2_id": "bin2_id"}
            )
            .loc[:, lambda df: ["bin1_id", "bin2_id", "count", *df.columns[3:]]]
            .sort_values(by=["bin1_id", "bin2_id"])
        )

        bins = (
            self.genomic_bins.df.rename(
                columns={"Chromosome": "chrom", "Start": "start", "End": "end"}
            )
            .sort_values("genomic_bin_id")
            .assign(bin_id=lambda df: df["genomic_bin_id"])
            .set_index("genomic_bin_id")
        )

        cooler.create_cooler(
            cooler_fn,
            bins=bins,
            pixels=pixels,
            metadata=metadata,
            mode="w" if not os.path.exists(store) else "a",
            columns=pixels.columns[2:],
            dtypes=dict(zip(pixels.columns[2:], ["float32"] * len(pixels.columns[2:]))),
            ensure_sorted=True,
            ordered=True,
        )

        return cooler_fn


def link_common_cooler_tables(clr: os.PathLike):
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


def get_merged_cooler_metadata(coolers: Iterable[os.PathLike]):
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


def merge_coolers(coolers: Tuple, output: os.PathLike):
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

        metadata = get_merged_cooler_metadata(cooler_uris)

        with h5py.File(output, mode="a") as dest:
            dest[viewpoint.replace("::", "/resolutions/")].attrs[
                "metadata"
            ] = ujson.encode(metadata)

    # Reduce space by linking common tables (bins, chroms)
    link_common_cooler_tables(output)
