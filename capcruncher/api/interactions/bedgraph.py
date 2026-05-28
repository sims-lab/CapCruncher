from __future__ import annotations

import re
from pathlib import Path
from types import NotImplementedType
from typing import Self

import cooler
import numpy as np
import pandas as pd
import pyranges1 as pr
from loguru import logger

from capcruncher.api.interactions.cooler.binning import CoolerBinner
from capcruncher.types import Normalisation
from capcruncher.utils import is_valid_bed

def _bedgraph_to_pyranges(bedgraph: pd.DataFrame) -> pr.PyRanges:
    return pr.PyRanges(
        bedgraph.rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
    )


def _pyranges_to_bedgraph(ranges: pr.PyRanges) -> pd.DataFrame:
    return ranges.rename(
        columns={"Chromosome": "chrom", "Start": "start", "End": "end"}
    )


def _cluster_multi_viewpoint_bins_bedgraph(bedgraph: pd.DataFrame) -> pd.DataFrame:
    gr_bdg = _bedgraph_to_pyranges(bedgraph)

    return (
        gr_bdg.cluster()
        .groupby("Cluster")
        .agg(
            {
                "count": "sum",
                "Start": "min",
                "End": "max",
                "Chromosome": "first",
            }
        )
        .reset_index()
        .rename(columns={"Start": "start", "End": "end", "Chromosome": "chrom"})[
            ["chrom", "start", "end", "count"]
        ]
    )


class CoolerBedGraph:
    """Generates a bedgraph file from a cooler file created by interactions-store.

    Attributes:
     cooler (cooler.Cooler): Cooler file to use for bedgraph production
     capture_name (str): Name of capture probe being processed.
     sparse (bool): Only output bins with interactions.
     only_cis (bool): Only output cis interactions.

    """

    def __init__(
        self,
        uri: str,
        sparse: bool = True,
        only_cis: bool = False,
        region_to_limit: str | None = None,
    ):
        """
        Args:
            uri (str): Path to cooler group in hdf5 file.
            sparse (bool, optional): Only output non-zero bins. Defaults to True.
        """
        self._sparse = sparse
        self._only_cis = only_cis

        logger.info(f"Reading {uri}")
        self._cooler = cooler.Cooler(uri)
        self.viewpoint_name = self._cooler.info["metadata"]["viewpoint_name"]
        self._viewpoint_bins = self._cooler.info["metadata"]["viewpoint_bins"]

        if len(self._viewpoint_bins) > 1:
            self.multiple_viewpoint_bins = True
            logger.warning(
                f"Viewpoint {self.viewpoint_name} has multiple bins! {self._viewpoint_bins}. Proceed with caution!"
            )
        else:
            self.multiple_viewpoint_bins = False

        self.viewpoint_chroms = self._cooler.info["metadata"]["viewpoint_chrom"]
        self.n_cis_interactions = self._cooler.info["metadata"]["n_cis_interactions"]
        logger.info(f"Processing {self.viewpoint_name}")

        if only_cis:
            pixels = []
            bins = []
            for chrom in self.viewpoint_chroms:
                _bins = self._cooler.bins().fetch(chrom)
                viewpoint_chrom_bins = self._bins["name"]
                _pixels = (
                    self._cooler.pixels()
                    .fetch(self.viewpoint_chroms)
                    .query(
                        "(bin1_id in @viewpoint_chrom_bins) and (bin2_id in @viewpoint_chrom_bins)"
                    )
                )
                _bins = self._cooler.bins().fetch(chrom)

                pixels.append(_pixels)
                bins.append(_bins)

            self._pixels = pd.concat(pixels)
            self._bins = pd.concat(bins)

        elif region_to_limit:
            self._pixels = self._cooler.pixels().fetch(region_to_limit)
            self._bins = self._cooler.bins().fetch(region_to_limit)

        else:
            self._pixels = self._cooler.pixels()[:]
            # TODO: Avoid this if possible as reading all bins into memory
            self._bins = self._cooler.bins()[:]

        # Ensure name column is present
        self._bins = (
            self._bins.assign(name=lambda df: df.index)
            if "name" not in self._bins.columns
            else self._bins
        )
        self._reporters = None

    def _get_reporters(self) -> pd.DataFrame:
        logger.info("Extracting reporters")
        concat_ids = pd.concat([self._pixels["bin1_id"], self._pixels["bin2_id"]])
        concat_ids_filt = concat_ids.loc[lambda ser: ser.isin(self._viewpoint_bins)]
        pixels = self._pixels.loc[concat_ids_filt.index]

        df_interactions = pd.DataFrame()
        df_interactions["capture"] = np.where(
            pixels["bin1_id"].isin(self._viewpoint_bins),
            pixels["bin1_id"],
            pixels["bin2_id"],
        )

        df_interactions["reporter"] = np.where(
            pixels["bin1_id"].isin(self._viewpoint_bins),
            pixels["bin2_id"],
            pixels["bin1_id"],
        )

        df_interactions["count"] = pixels["count"].values

        return df_interactions.sort_values(["capture", "reporter"]).reset_index(
            drop=True
        )

    def extract_bedgraph(
        self, normalisation: Normalisation = Normalisation.RAW, **norm_kwargs
    ) -> pd.DataFrame:
        logger.info("Generating bedgraph")
        df_bdg = (
            self._bins.merge(
                self.reporters,
                left_on="name",
                right_on="reporter",
                how="inner" if self._sparse else "outer",
            )[["chrom", "start", "end", "count"]]
            .assign(count=lambda df: df["count"].fillna(0))
            .sort_values(["chrom", "start"])
        )

        if self.multiple_viewpoint_bins:
            df_bdg = _cluster_multi_viewpoint_bins_bedgraph(df_bdg)

        if normalisation != Normalisation.RAW:
            logger.info("Normalising bedgraph")
            self._normalise_bedgraph(df_bdg, method=normalisation, **norm_kwargs)

        return df_bdg

    @property
    def reporters(self) -> pd.DataFrame:
        """Interactions with capture fragments/bins.

        Returns:
         pd.DataFrame: DataFrame containing just bins interacting with the capture probe.
        """

        if self._reporters is not None:
            return self._reporters
        else:
            self._reporters = self._get_reporters()
            return self._reporters

    def _normalise_bedgraph(
        self,
        bedgraph: pd.DataFrame,
        scale_factor: float = 1e6,
        method: Normalisation = Normalisation.N_CIS,
        region: Path | str | None = None,
    ) -> None:
        """Normalises the bedgraph (in place).

        Uses the number of cis interactions to normalise the bedgraph counts.

        Args:
         scale_factor (int, optional): Scaling factor for normalisation. Defaults to 1e6.

        Returns:
         pd.DataFrame: Normalised bedgraph formatted DataFrame
        """

        if method == Normalisation.RAW:
            pass
        elif method == Normalisation.N_CIS:
            self._normalise_by_n_cis(bedgraph, scale_factor)
        elif method == Normalisation.REGION:
            if region is None:
                raise ValueError("Region based normalisation requires a BED file.")
            self._normalise_by_regions(bedgraph, scale_factor, region)

    def _normalise_by_n_cis(
        self, bedgraph: pd.DataFrame, scale_factor: float
    ) -> None:
        bedgraph["count"] = (bedgraph["count"] / self.n_cis_interactions) * scale_factor

    def _normalise_by_regions(
        self, bedgraph: pd.DataFrame, scale_factor: float, regions: Path | str
    ) -> None:
        if not is_valid_bed(regions):
            raise ValueError(
                "A valid bed file is required for region based normalisation"
            )

        df_viewpoint_norm_regions = pd.read_csv(
            regions, sep="\t", names=["chrom", "start", "end", "name"]
        )
        df_viewpoint_norm_regions = df_viewpoint_norm_regions.loc[
            lambda df: df["name"].str.contains(self.viewpoint_name)
        ]

        gr_bedgraph = _bedgraph_to_pyranges(
            bedgraph.reset_index(names="bedgraph_row_id")
        )
        gr_regions = _bedgraph_to_pyranges(
            df_viewpoint_norm_regions.rename(columns={"name": "region_name"})
        )
        df_counts_in_regions = gr_bedgraph.join_overlaps(
            gr_regions, strand_behavior="ignore"
        ).drop_duplicates("bedgraph_row_id")
        total_counts_in_region = df_counts_in_regions["count"].sum()

        bedgraph["count"] = (bedgraph["count"] / total_counts_in_region) * scale_factor

    def to_pyranges(
        self, normalisation: Normalisation = Normalisation.RAW, **norm_kwargs
    ) -> pr.PyRanges:
        return pr.PyRanges(
            self.extract_bedgraph(
                normalisation=normalisation, **norm_kwargs
            ).rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
        )


class CoolerBedGraphWindowed(CoolerBedGraph):
    def __init__(
        self,
        cooler_fn: str,
        binsize: int = 5_000,
        binner: CoolerBinner | None = None,
        sparse: bool = True,
    ):
        super().__init__(cooler_fn, sparse=sparse)

        self.cooler = cooler.Cooler(cooler_fn)
        self.binner = binner if binner else CoolerBinner(cooler_fn, binsize=binsize)
        self.binsize = self.binner.binsize
        self._bins_genomic = self.binner.bins
        self.capture_bins = self.cooler.info["metadata"]["capture_bins"]
        self.capture_bins_genomic = self.binner.viewpoint_bins

    def _get_bedgraph(self):
        bct = self.binner.bin_conversion_table
        reporters = self.reporters

        # Merge reporters with windows
        bedgraph_frag = bct.merge(
            reporters, left_on="name_fragment", right_on="reporter"
        ).drop(columns=["capture", "name_fragment"])

        # Get aggregated count
        count_aggregated = (
            bedgraph_frag.groupby("name_bin").agg({"count": "sum"}).reset_index()
        )

        # Merge bins with aggregated counts
        bedgraph_bins = self.binner.bins.merge(
            count_aggregated,
            left_on="name",
            right_on="name_bin",
            how="inner" if self._sparse else "outer",
        ).drop(columns=["name_bin"])[["chrom", "start", "end", "count"]]

        return bedgraph_bins

    def _normalise_bedgraph(
        self,
        bedgraph: pd.DataFrame,
        scale_factor: float = 1e6,
        method: Normalisation = Normalisation.N_CIS,
        region: Path | str | None = None,
    ) -> None:
        bct = self.binner.bin_conversion_table
        reporters = self.reporters

        bedgraph_frag = (
            bct.merge(reporters, left_on="name_fragment", right_on="reporter")
            .drop(columns=["capture", "name_fragment"])
            .assign(
                count_overfrac_norm=lambda df: df["count"] * df["overlap_fraction"],
                count_overfrac_n_interact_norm=lambda df: (
                    df["count_overfrac_norm"] / self.n_cis_interactions
                )
                * scale_factor,
            )
        )

        count_aggregated = (
            bedgraph_frag.groupby("name_bin")
            .agg({"count_overfrac_n_interact_norm": "mean"})
            .reset_index()
        )

        bedgraph_bins = self.binner.bins.merge(
            count_aggregated,
            left_on="name",
            right_on="name_bin",
            how="inner" if self._sparse else "outer",
        ).drop(columns=["name_bin"])[
            ["chrom", "start", "end", "count_overfrac_n_interact_norm"]
        ]

        bedgraph_bins.columns = ["chrom", "start", "end", "count"]

        bedgraph[["chrom", "start", "end", "count"]] = bedgraph_bins[
            ["chrom", "start", "end", "count"]
        ]

    @property
    def reporters_binned(self):
        reporters = self.reporters
        reporters_binned = (
            reporters.merge(
                self.binner.bin_conversion_table[["name_bin", "name_fragment"]],
                left_on="capture",
                right_on="name_fragment",
            )
            .merge(
                self.binner.bin_conversion_table[["name_bin", "name_fragment"]],
                left_on="reporter",
                right_on="name_fragment",
            )
            .groupby(["name_bin_x", "name_bin_y"])
            .agg({"count": "sum"})
            .reset_index()[["name_bin_x", "name_bin_y", "count"]]
        )
        reporters_binned.columns = ["capture", "reporter", "count"]
        return reporters_binned


class CCBedgraph:
    def __init__(
        self,
        path: Path | str | None = None,
        df: pd.DataFrame | None = None,
        capture_name: str = "",
        capture_chrom: str = "",
        capture_start: str = "",
        capture_end: str = "",
    ) -> None:
        self.fn = path
        self.df = df

        if self.fn:
            self.df = pd.read_csv(
                self.fn, sep="\t", header=None, names=["chrom", "start", "end", "score"]
            )

        self.capture_name = capture_name
        self.capture_chrom = capture_chrom
        self.capture_start = capture_start
        self.capture_end = capture_end

    @property
    def _df(self) -> pd.DataFrame:
        if self.df is None:
            raise ValueError("CCBedgraph requires either a path or dataframe.")
        return self.df

    @property
    def score(self) -> pd.Series:
        return self._df.rename(columns={"score": self.fn})[self.fn]

    @property
    def coordinates(self) -> pd.DataFrame:
        return self._df.loc[:, "chrom":"end"]

    def to_pyranges(self) -> pr.PyRanges:
        return self._df.rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        ).pipe(pr.PyRanges)

    def to_file(self, path: Path | str) -> None:
        self._df.to_csv(path, sep="\t", header=None, index=None)

    def __add__(self, other: object) -> Self | NotImplementedType:
        if isinstance(other, CCBedgraph):
            self._df["score"] = self._df["score"] + other._df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self._df["score"] = self._df["score"] + other
            return self

        else:
            return NotImplemented

    def __sub__(self, other: object) -> Self | NotImplementedType:
        if isinstance(other, CCBedgraph):
            self._df["score"] = self._df["score"] - other._df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self._df["score"] = self._df["score"] - other
            return self

        else:
            return NotImplemented

    def __mul__(self, other: object) -> Self | NotImplementedType:
        if isinstance(other, CCBedgraph):
            self._df["score"] = self._df["score"] * other._df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self._df["score"] = self._df["score"] * other
            return self

        else:
            return NotImplemented

    def __truediv__(self, other: object) -> Self | NotImplementedType:
        if isinstance(other, CCBedgraph):
            self._df["score"] = self._df["score"] / other._df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self._df["score"] = self._df["score"] / other
            return self

        else:
            return NotImplemented


def cooler_to_bedgraph(
    clr: str,
    regions_of_interest: Path | str | None = None,
    viewpoint_distance: int | None = None,
    **kwargs,
) -> pd.DataFrame:
    if viewpoint_distance:
        viewpoint_coords = cooler.Cooler(clr).info["metadata"]["viewpoint_coords"][0]
        viewpoint_coords = re.split("[:-]", viewpoint_coords)

        viewpoint_chrom = viewpoint_coords[0]
        viewpoint_start = max(0, int(viewpoint_coords[1]) - viewpoint_distance)
        viewpoint_chromsize = cooler.Cooler(clr).chromsizes[viewpoint_chrom]
        viewpoint_end = min(
            int(viewpoint_coords[1]) + viewpoint_distance, viewpoint_chromsize
        )
        region_to_limit = f"{viewpoint_chrom}:{viewpoint_start}-{viewpoint_end}"
        logger.info(
            f"Limiting bedgraph to {viewpoint_chrom}:{viewpoint_start}-{viewpoint_end}"
        )
    else:
        region_to_limit = None

    norm = kwargs.get("normalisation", "raw")
    scale_factor = kwargs.get("scale_factor", 1e6)
    region = kwargs.get("region", None)

    bedgraph = CoolerBedGraph(clr, region_to_limit=region_to_limit).extract_bedgraph(
        normalisation=norm,
        scale_factor=scale_factor,
        region=region,
    )

    if regions_of_interest:
        pr_roi = pr.read_bed(Path(regions_of_interest))
        pr_bedgraph = _bedgraph_to_pyranges(bedgraph)
        pr_bedgraph = pr_bedgraph.join_overlaps(pr_roi, strand_behavior="same")

        bedgraph = _pyranges_to_bedgraph(pr_bedgraph)[
            ["chrom", "start", "end", "count"]
        ]

    return bedgraph
