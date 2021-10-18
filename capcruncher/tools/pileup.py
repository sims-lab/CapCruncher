import pandas as pd
import numpy as np
from pybedtools import BedTool
import cooler
from typing import Union
from capcruncher.tools.storage import CoolerBinner
from capcruncher.utils import is_valid_bed
import os


class CoolerBedGraph:
    """Generates a bedgraph file from a cooler file created by interactions-store.

    Attributes:
     cooler (cooler.Cooler): Cooler file to use for bedgraph production
     capture_name (str): Name of capture probe being processed.
     sparse (bool): Only output bins with interactions.
     only_cis (bool): Only output cis interactions.

    """

    def __init__(self, cooler_fn: str, sparse: bool = True, only_cis: bool = False):
        """
        Args:
            cooler_fn (str): Path to cooler group in hdf5 file.
            sparse (bool, optional): Only output non-zero bins. Defaults to True.
        """
        self.sparse = sparse
        self.only_cis = only_cis

        self.cooler = cooler.Cooler(cooler_fn)
        self.viewpoint_name = self.cooler.info["metadata"]["viewpoint_name"]
        self._viewpoint_bins = self.cooler.info["metadata"]["viewpoint_bins"]
        self._viewpoint_chrom = self.cooler.info["metadata"]["viewpoint_chrom"]
        self._n_cis_interactions = self.cooler.info["metadata"]["n_cis_interactions"]
        self._bins = self.cooler.bins()[:]
        self._pixels = (
            self.cooler.pixels()[:]
            if not only_cis
            else self.cooler.pixels().fetch(self._viewpoint_chrom)
        )

        self._bedgraph = None
        self._reporters = None

    def _get_reporters(self):

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

    def _get_bedgraph(self):

        df_bdg = (
            self._bins.merge(
                self.reporters,
                left_on="name",
                right_on="reporter",
                how="inner" if self.sparse else "outer",
            )[["chrom", "start", "end", "count"]]
            .assign(count=lambda df: df["count"].fillna(0))
            .sort_values(["chrom", "start"])
        )

        return df_bdg

    @property
    def bedgraph(self) -> pd.DataFrame:
        """
        Returns:
         pd.DataFrame: DataFrame in bedgraph format.
        """
        if self._bedgraph is not None:
            return self._bedgraph
        else:
            self._bedgraph = self._get_bedgraph()
            return self._bedgraph

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

    def normalise_bedgraph(
        self, scale_factor=1e6, method: str = "n_cis", region: str = None
    ):
        """Normalises the bedgraph.

        Uses the number of cis interactions to normalise the bedgraph counts.

        Args:
         scale_factor (int, optional): Scaling factor for normalisation. Defaults to 1e6.

        Returns:
         pd.DataFrame: Normalised bedgraph formatted DataFrame
        """

        if method == "raw":
            pass
        elif method == "n_cis":
            self._normalise_by_n_cis(scale_factor)
        elif method == "region":
            self._normalise_by_regions(scale_factor, region)

    def _normalise_by_n_cis(self, scale_factor: float):
        self.bedgraph["count"] = (
            self.bedgraph["count"] / self._n_cis_interactions
        ) * scale_factor

    def _normalise_by_regions(self, scale_factor: float, regions: str):

        if not is_valid_bed(regions):
            raise ValueError(
                "A valid bed file is required for region based normalisation"
            )

        df_viewpoint_norm_regions = pd.read_csv(regions, sep="\t", names=["chrom", "start", "end", "name"])
        df_viewpoint_norm_regions = df_viewpoint_norm_regions.loc[lambda df: df["name"].str.contains(self.viewpoint_name)]

        counts_in_regions = []
        for region in df_viewpoint_norm_regions.itertuples():
            counts_in_regions.append(
                self.bedgraph.query(
                    "(chrom == @region.chrom) and (start >= @region.start) and (start <= @region.end)"
                )
            )
        
        df_counts_in_regions = pd.concat(counts_in_regions)
        total_counts_in_region = df_counts_in_regions['count'].sum()

        self.bedgraph["count"] = (
            self.bedgraph["count"] / total_counts_in_region
        ) * scale_factor
    

    def to_file(self, fn: os.PathLike):
        """Outputs the bedgraph dataframe to a file.

        If normalise is True, will also normalise the counts by the number of cis interactions.

        Args:
         fn (os.PathLike): Output file name.
         normalise (bool, optional): Normalise the bedgraph before writing to file. Defaults to False.
        """
        self.bedgraph.to_csv(fn, sep="\t", header=None, index=False)


class CoolerBedGraphWindowed(CoolerBedGraph):
    def __init__(
        self,
        cooler_fn: str,
        binsize: int = 5e3,
        binner: CoolerBinner = None,
        sparse=True,
    ):

        super(CoolerBedGraphWindowed, self).__init__(cooler_fn, sparse=sparse)

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
            how="inner" if self.sparse else "outer",
        ).drop(columns=["name_bin"])[["chrom", "start", "end", "count"]]

        return bedgraph_bins

    def normalise_bedgraph(self, scale_factor=1e6):

        bct = self.binner.bin_conversion_table
        reporters = self.reporters

        bedgraph_frag = (
            bct.merge(reporters, left_on="name_fragment", right_on="reporter")
            .drop(columns=["capture", "name_fragment"])
            .assign(
                count_overfrac_norm=lambda df: df["count"] * df["overlap_fraction"],
                count_overfrac_n_interact_norm=lambda df: (
                    df["count_overfrac_norm"] / self._n_cis_interactions
                )
                * scale_factor,
            ),
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
            how="inner" if self.sparse else "outer",
        ).drop(columns=["name_bin"])[
            ["chrom", "start", "end", "count_overfrac_n_interact_norm"]
        ]

        bedgraph_bins.columns = ["chrom", "start", "end", "count"]

        return bedgraph_bins

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


class CCBedgraph(object):
    def __init__(
        self,
        path=None,
        df=None,
        capture_name="",
        capture_chrom="",
        capture_start="",
        capture_end="",
    ):

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
    def score(self):
        return self.df.rename(columns={"score": self.fn})[self.fn]

    @property
    def coordinates(self):
        return self.df.loc[:, "chrom":"end"]

    def to_bedtool(self):
        return self.df.pipe(BedTool.from_dataframe)

    def to_file(self, path):
        self.df.to_csv(path, sep="\t", header=None, index=None)

    def __add__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] + other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] + other
            return self

        else:
            return NotImplementedError()

    def __sub__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] - other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] - other
            return self

        else:
            return NotImplementedError()

    def __mul__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] * other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] * other
            return self

        else:
            return NotImplementedError()

    def __truediv__(self, other):
        if isinstance(other, CCBedgraph):
            self.df["score"] = self.df["score"] / other.df["score"]
            return self

        elif isinstance(other, (np.ndarray, pd.Series, int, float)):
            self.df["score"] = self.df["score"] / other
            return self

        else:
            return NotImplementedError()
