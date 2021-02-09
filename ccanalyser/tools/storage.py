import os
import sys
import pandas as pd
import numpy as np
from pybedtools import BedTool
import cooler
import h5py
from joblib import Parallel, delayed
from typing import Union
from ccanalyser.utils import split_intervals_on_chrom, intersect_bins

# Required for initial storage


def get_capture_coords(oligo_file: str, oligo_name: str):
    df_oligos = BedTool(oligo_file).to_dataframe()
    df_oligos = df_oligos.query(f'name == "oligo_name"')
    return BedTool.from_dataframe(df_oligos)


def get_capture_bins_from_oligos(
    oligo_name: str, oligo_file: str, fragments: Union[str, pd.DataFrame]
):

    bt_oligo = get_capture_coords(oligo_file, oligo_name)

    if isinstance(fragments, str):
        bt_rf = BedTool(fragments)
    elif isinstance(fragments, pd.DataFrame):
        bt_rf = BedTool.from_dataframe(fragments)
    else:
        raise ValueError("Provide fragments as either a filename string or DataFrame")

    bt_bins = bt_rf.intersect(bt_oligo, f=0.51, sorted=True)
    return bt_bins.to_dataframe()["name"].to_list()


def create_cooler_cc(
    fn: str,
    bins: pd.DataFrame,
    pixels: pd.DataFrame,
    capture_name: str,
    capture_oligos: str,
    capture_bins: Union[int, list] = None,
    **cooler_kwargs,
):

    capture_coords = "\n".join(
        ["\t".join(x for x in get_capture_coords(capture_oligos, capture_name))]
    )

    if not capture_bins:
        capture_bins = get_capture_bins_from_oligos(capture_name, capture_oligos, bins)

    metadata = {
        "capture_bins": [int(x) for x in (capture_bins,)],
        "capture_name": capture_name,
        "capture_coords": capture_coords,
    }

    # Create cooler
    cooler.create_cooler(
        f"{fn}::{capture_name}",
        bins=bins,
        pixels=pixels,
        metadata=metadata,
        mode="w" if not os.path.exists(fn) else "a",
        **cooler_kwargs,
    )

    return f"{fn}/{capture_name}"


# Required for binning interactions into genomic regions


class CoolerBinner:
    def __init__(
        self, cooler_fn, binsize=5e4, normalise=True, scale_factor=1e6, n_cores=8
    ):

        self.cooler = cooler.Cooler(cooler_fn)
        self.binsize = binsize
        self.normalise = normalise
        self.scale_factor = scale_factor
        self.n_cores = n_cores

        self.bins_fragments = self.cooler.bins()[:]
        self.bins_genomic = self._get_bins()
        self._bin_conversion_table = None
        self._pixel_conversion_table = None
        self._pixels = None

    def _get_bins(self):
        return (
            cooler.util.binnify(chromsizes=self.cooler.chromsizes, binsize=self.binsize)
            .reset_index()
            .rename(columns={"index": "name"})[["chrom", "start", "end", "name"]]
            .assign(
                start=lambda df: df["start"].astype(int),
                end=lambda df: df["end"].astype(int),
            )
        )

    def _get_bin_conversion_table(self):

        bins_genomic_by_chrom = split_intervals_on_chrom(self.bins_genomic)
        bins_fragments_by_chrom = split_intervals_on_chrom(self.bins_fragments)

        bins_intersections = Parallel(n_jobs=self.n_cores)(
            delayed(intersect_bins)(
                bins_genomic_by_chrom[chrom], bins_fragments_by_chrom[chrom]
            )
            for chrom in bins_genomic_by_chrom
        )

        df_bins_intersections = pd.concat(bins_intersections, ignore_index=True).rename(
            columns=lambda c: c.replace("_1", "_bin").replace("_2", "_fragment")
        )

        df_bins_intersections["overlap_fraction"] = df_bins_intersections["overlap"] / (
            df_bins_intersections["end_fragment"]
            - df_bins_intersections["start_fragment"]
        )

        return df_bins_intersections

    def _get_pixel_conversion_table(self):

        pixels = self.cooler.pixels()[:]
        pixels_conv = (
            pixels.merge(
                self.bin_conversion_table[["name_bin", "name_fragment"]].add_suffix(
                    "_1"
                ),
                left_on="bin1_id",
                right_on="name_fragment_1",
            )
            .merge(
                self.bin_conversion_table[["name_bin", "name_fragment"]].add_suffix(
                    "_2"
                ),
                left_on="bin2_id",
                right_on="name_fragment_2",
            )
            .drop(columns=["name_fragment_1", "name_fragment_2"])
        )

        return pixels_conv

    def _get_pixels(self):

        df_pixels = (
            self.pixel_conversion_table.groupby(["name_bin_1", "name_bin_2"])["count"]
            .sum()
            .to_frame()
            .reset_index()
        )

        df_pixels.columns = ["bin1_id", "bin2_id", "count"]
        df_pixels = df_pixels.loc[lambda df: df["bin1_id"] != df["bin2_id"]]

        return df_pixels

    @property
    def bins(self):
        return self.bins_genomic

    @property
    def bin_conversion_table(self):
        if self._bin_conversion_table is not None:
            return self._bin_conversion_table
        else:
            self._bin_conversion_table = self._get_bin_conversion_table()
            return self._bin_conversion_table

    @property
    def capture_bins(self):
        capture_frags = self.cooler.info["metadata"]["capture_bins"]
        return self.bin_conversion_table.loc[
            lambda df: df["name_fragment"].isin(capture_frags)
        ]["name_bin"].values

    @property
    def pixel_conversion_table(self):
        if self._pixel_conversion_table is not None:
            return self._pixel_conversion_table
        else:
            self._pixel_conversion_table = self._get_pixel_conversion_table()
            return self._pixel_conversion_table

    @property
    def pixels(self):
        if self._pixels is not None:
            return self._pixels
        else:
            self._pixels = self._get_pixels()
            return self._pixels

    def normalise_pixels(
        self,
        n_fragment_correction=True,
        n_interaction_correction=True,
        scale_factor=1e6,
    ):

        total_interactions = self.pixels["count"].sum()

        if n_fragment_correction:
            df_nrf = (
                self.bin_conversion_table.groupby("name_bin").size().to_frame("n_rf")
            )
            pixels = self.pixels.merge(
                df_nrf.add_prefix("bin1_"), left_on="bin1_id", right_index=True
            ).merge(df_nrf.add_prefix("bin2_"), left_on="bin2_id", right_index=True)

            self.pixels["count_n_rf_norm"] = self.pixels["count"] / (
                pixels["bin1_n_rf"] * pixels["bin2_n_rf"]
            )

        if n_interaction_correction:
            self.pixels["count_n_interactions_norm"] = (
                self.pixels["count"] / total_interactions
            ) * scale_factor

        if n_fragment_correction and n_interaction_correction:

            self.pixels["count_n_rf_n_interactions_norm"] = (
                self.pixels["count_n_rf_norm"] / total_interactions
            ) * scale_factor

    def to_cooler(self, store, normalise=False, **normalise_options):

        capture_bins = self.capture_bins
        capture_name = self.cooler.info["metadata"]["capture_name"]
        capture_coords = self.cooler.info['metadata']['capture_coords']

        metadata = {
            "capture_bins": [int(x) for x in (capture_bins,)],
            "capture_name": capture_name,
            "capture_coords":  capture_coords
        }

        if normalise:
            self.normalise_pixels(**normalise_options)

        # Create cooler
        cooler.create_cooler(
            f"{store}::{capture_name}/resolutions/{int(self.binsize)}",
            bins=self.bins,
            pixels=self.pixels,
            metadata=metadata,
            mode="w" if not os.path.exists(store) else "a",
            columns=self.pixels.columns[2:],
        )

        return f"{store}/{capture_name}/resolutions/{int(self.binsize)}"
