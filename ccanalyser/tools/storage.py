import os
import sys
import pandas as pd
import numpy as np
from pybedtools import BedTool
import cooler
import h5py
from joblib import Parallel, delayed
from typing import Union
from natsort import natsorted
from ccanalyser.utils import split_intervals_on_chrom, intersect_bins

# Required for initial storage


def get_capture_coords(oligo_file: str, oligo_name: str):
    df_oligos = BedTool(oligo_file).to_dataframe()
    df_oligos = df_oligos.query(f'name == "{oligo_name}"')
    return df_oligos.iloc[0]


# def get_capture_bins_from_oligos(
#     oligo_name: str, oligo_file: str, fragments: Union[str, pd.DataFrame]
# ):

#     bt_oligo = get_capture_coords(oligo_file, oligo_name)

#     if isinstance(fragments, str):
#         bt_rf = BedTool(fragments)
#     elif isinstance(fragments, pd.DataFrame):
#         bt_rf = BedTool.from_dataframe(fragments)
#     else:
#         raise ValueError("Provide fragments as either a filename string or DataFrame")

#     bt_bins = bt_rf.intersect(bt_oligo, f=0.51, sorted=True)
#     return bt_bins.to_dataframe()["name"].to_list()


def get_capture_bins(bins, oligo_chrom, oligo_start, oligo_end):

    return bins.query(
        f'chrom == "{oligo_chrom}" and start >= {oligo_start} and end <= {oligo_end}'
    )["name"]


def create_cooler_cc(
    output_prefix: str,
    bins: pd.DataFrame,
    pixels: pd.DataFrame,
    capture_name: str,
    capture_oligos: str,
    capture_bins: Union[int, list] = None,
    suffix=None,
    **cooler_kwargs,
):

    capture_coords = get_capture_coords(capture_oligos, capture_name)

    if capture_coords is None:
        raise ValueError(f"Incorrect capture name specified: {capture_name}.")

    if not capture_bins:
        capture_bins = get_capture_bins(
            bins,
            capture_coords["chrom"],
            capture_coords["start"],
            capture_coords["end"],
        )
        capture_bins = [int(x) for x in capture_bins]

    elif isinstance(capture_bins, int):
        capture_bins = [
            int(capture_bins),
        ]

    elif isinstance(capture_bins, (np.array, pd.Series)):
        capture_bins = [int(x) for x in capture_bins]

    metadata = {
        "capture_bins": capture_bins,
        "capture_name": capture_name,
        "capture_coords": f'{capture_coords["chrom"]}:{capture_coords["start"]}-{capture_coords["end"]}',
    }


    if os.path.exists(output_prefix): # Will append to a prexisting file if one is supplied
        cooler_fn = f"{output_prefix}::/{capture_name}"
    else:
        cooler_fn = f"{output_prefix.replace('.hdf5', '')}.{capture_name}.{suffix if suffix else ''}.hdf5"
        
    
    cooler.create_cooler(
        cooler_fn,
        bins=bins,
        pixels=pixels,
        metadata=metadata,
        mode="w" if not os.path.exists(cooler_fn.split('::')[0]) else "a",
        **cooler_kwargs,
    )

    return cooler_fn


class GenomicBinner:
    def __init__(self, chromsizes, fragments, binsize=5000, n_cores=8, overlap_fraction=0.2):

        self.chromsizes = self._format_chromsizes(chromsizes)
        self.fragments = self._format_fragments(fragments)
        self.binsize = binsize
        self.overlap_fraction = overlap_fraction

        self.bins_genomic = self._get_bins()
        self._bin_conversion_table = None
        self.n_cores = n_cores

    def _get_bins(self):
        return (
            cooler.util.make_bintable(chromsizes=self.chromsizes, binsize=self.binsize)
            .reset_index()
            .rename(columns={"index": "name"})[["chrom", "start", "end", "name"]]
            .assign(
                start=lambda df: df["start"].astype(int),
                end=lambda df: df["end"].astype(int),
            )
        )

    def _format_chromsizes(self, chromsizes):

        _chromsizes = pd.Series()

        if isinstance(chromsizes, str):
            _chromsizes = pd.read_csv(
                chromsizes, sep="\t", index_col=0, header=None
            ).loc[lambda df: natsorted(df.index), 0]

        elif isinstance(chromsizes, pd.DataFrame):
            if chromsizes.index.astype(str).str.contains("^chr.*"):
                _chromsizes = chromsizes.iloc[:, 0]

        elif isinstance(chromsizes, pd.Series):
            _chromsizes = chromsizes

        if not _chromsizes.empty:
            return _chromsizes
        else:
            raise ValueError("Chromsizes supplied in the wrong format")

    def _natsort_dataframe(self, df, column):

        df_by_key = {k: df for k, df in df.groupby(column)}

        _df = pd.DataFrame()
        for k in natsorted(df_by_key):
            if _df is not None:
                _df = pd.concat([_df, df_by_key[k]])
            else:
                _df = df_by_key[k]

        return _df

    def _format_fragments(self, fragments):

        if isinstance(fragments, str):
            _fragments = pd.read_csv(
                fragments,
                sep="\t",
                index_col=0,
                header=None,
                names=["chrom", "start", "end", "name"],
            )
        elif isinstance(fragments, pd.DataFrame):
            _fragments = fragments

        return self._natsort_dataframe(_fragments, "chrom")

    def _get_bin_conversion_table(self):

        bins_genomic_by_chrom = split_intervals_on_chrom(self.bins_genomic)
        bins_fragments_by_chrom = split_intervals_on_chrom(self.fragments)

        shared_chroms = set(bins_fragments_by_chrom) & set(bins_genomic_by_chrom)

        bins_intersections = Parallel(n_jobs=self.n_cores)(
            delayed(intersect_bins)(
                bins_genomic_by_chrom[chrom], bins_fragments_by_chrom[chrom], F=self.overlap_fraction,
            )
            for chrom in natsorted(shared_chroms)
        )

        df_bins_intersections = pd.concat(bins_intersections, ignore_index=True).rename(
            columns=lambda c: c.replace("_1", "_bin").replace("_2", "_fragment")
        )

        df_bins_intersections["overlap_fraction"] = df_bins_intersections["overlap"] / (
            df_bins_intersections["end_fragment"]
            - df_bins_intersections["start_fragment"]
        )

        return df_bins_intersections

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


class CoolerBinner:
    def __init__(
        self,
        cooler_fn,
        binsize=None,
        scale_factor=1e6,
        n_cores=8,
        binner=None,
    ):

        self.cooler = cooler.Cooler(cooler_fn)
        self.binner = binner or GenomicBinner(
            chromsizes=self.cooler.chromsizes,
            fragments=self.cooler.bins()[:],
            binsize=binsize,
        )

        self.binsize = self.binner.binsize
        self.scale_factor = scale_factor
        self.n_cores = n_cores

        self.bins_fragments = self.cooler.bins()[:]

        self._bin_conversion_table = None
        self._pixel_conversion_table = None
        self._pixels = None

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
        df_pixels["bin1_id_corrected"] = np.where(
            df_pixels["bin1_id"] > df_pixels["bin2_id"],
            df_pixels["bin2_id"],
            df_pixels["bin1_id"],
        )
        df_pixels["bin2_id_corrected"] = np.where(
            df_pixels["bin1_id"] > df_pixels["bin2_id"],
            df_pixels["bin1_id"],
            df_pixels["bin2_id"],
        )
        df_pixels = df_pixels.loc[lambda df: df["bin1_id"] != df["bin2_id"]]
        df_pixels = df_pixels.loc[:, 
            ["bin1_id_corrected", "bin2_id_corrected", "count"]
        ].rename(columns=lambda col: col.replace("_corrected", ""))

        return df_pixels

    @property
    def bins(self):
        return self.binner.bins

    @property
    def bin_conversion_table(self):
        return self.binner.bin_conversion_table

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
                self.pixels["count"] / scale_factor
            ) * total_interactions

        if n_fragment_correction and n_interaction_correction:

            self.pixels["count_n_rf_n_interactions_norm"] = (
                self.pixels["count_n_rf_norm"] / scale_factor
            ) * total_interactions

    def to_cooler(self, store, normalise=False, **normalise_options):

        capture_bins = self.capture_bins
        capture_name = self.cooler.info["metadata"]["capture_name"]
        capture_coords = self.cooler.info["metadata"]["capture_coords"]

        metadata = {
            "capture_bins": [int(x) for x in self.capture_bins],
            "capture_name": capture_name,
            "capture_coords": capture_coords,
        }

        if normalise:
            self.normalise_pixels(**normalise_options)

        if os.path.exists(store): # Will append to a prexisting file if one is supplied
            cooler_fn = f"{store}::/{capture_name}/resolutions/{self.binsize}"
        else:
            cooler_fn = f"{store.replace('.hdf5', '')}.{capture_name}.{self.binsize}.hdf5"
        
    
        cooler.create_cooler(
            cooler_fn,
            bins=self.bins,
            pixels=self.pixels,
            metadata=metadata,
            mode="w" if not os.path.exists(store) else "a",
            columns=self.pixels.columns[2:],
        )

        return cooler_fn
        
        
        
        
