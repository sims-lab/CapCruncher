import pandas as pd
import numpy as np
from pybedtools import BedTool
import cooler
from typing import Union


class CoolerBedGraph:
    def __init__(self, cooler_fn: str):

        self.cooler = cooler.Cooler(cooler_fn)
        self.bins = self.cooler.bins()[:]
        self.pixels = self.cooler.pixels()[:]
        self.capture_name = self.ccooler.info["metadata"]["capture_name"]
        self.capture_bins = self.cooler.info["metadata"]["capture_bins"]
        self._bedgraph = None
        self._reporters = None

    def _get_reporters(self):
        concat_ids = pd.concat([self.pixels["bin1_id"], self.pixels["bin2_id"]])
        concat_ids_filt = concat_ids.loc[lambda ser: ser.isin(self.capture_bins)]
        pixels = self.pixels.loc[concat_ids_filt.index]

        df_interactions = pd.DataFrame()
        df_interactions["capture"] = np.where(
            pixels["bin1_id"].isin(self.capture_bins),
            pixels["bin1_id"],
            pixels["bin2_id"],
        )

        df_interactions["reporter"] = np.where(
            pixels["bin1_id"].isin(self.capture_bins),
            pixels["bin2_id"],
            pixels["bin1_id"],
        )

        df_interactions["count"] = pixels["count"].values

        return df_interactions.sort_values(["capture", "reporter"]).reset_index(
            drop=True
        )

    def _get_bedgraph(self):

        df_bdg = (
            self.bins.merge(
                self.reporters, left_on="name", right_on="reporter", how="outer"
            )[["chrom", "start", "end", "count"]]
            .assign(count=lambda df: df["count"].fillna(0))
            .sort_values(["chrom", "start"])
        )

        return df_bdg

    @property
    def capture_chrom(self):
        return self.bins.loc[lambda df: df["name"].isin(self.capture_bins)][
            "chrom"
        ].mode()[0]

    @property
    def n_cis_interactions(self):
        return self.bedgraph.query(f'chrom == "{self.capture_chrom}"')["count"].sum()

    @property
    def bedgraph(self):
        if self._bedgraph is not None:
            return self._bedgraph
        else:
            self._bedgraph = self._get_bedgraph()
            return self._bedgraph

    @property
    def reporters(self):
        if self._reporters is not None:
            return self._reporters
        else:
            self._reporters = self._get_reporters()
            return self._reporters

    def normalise_bedgraph(self, scale_factor=1e6):
        df_bdg = self.bedgraph
        df_bdg["count"] = (df_bdg["count"] / self.n_cis_interactions) * scale_factor
        return df_bdg

    def to_file(self, fn, normalise=False, **normalise_kwargs):

        if not normalise:
            self.bedgraph.to_csv(fn, sep="\t", header=None, index=False)
        else:
            self.normalise_bedgraph(**normalise_kwargs).to_csv(
                fn, sep="\t", header=None, index=False
            )


class CoolerBedGraphWindowed(CoolerBedGraph):
    def __init__(self, cooler_fn: str, binsize=5e3, binner=None):

        super(CoolerBedGraphWindowed, self).__init__(cooler_fn)

        self.cooler = cooler.Cooler(cooler_fn)
        self.binner = binner if binner else CoolerBinner(cooler_fn, binsize=binsize)
        self.binsize = self.binner.binsize
        self._bins_genomic = self.binner.bins
        self.capture_bins = self.cooler.info["metadata"]["capture_bins"]
        self.capture_bins_genomic = self.binner.capture_bins

    def _get_bedgraph(self):

        bct = self.binner.bin_conversion_table
        reporters = self.reporters

        bedgraph_frag = bct.merge(
            reporters, left_on="name_fragment", right_on="reporter"
        ).drop(columns=["capture", "name_fragment"])

        count_aggregated = (
            bedgraph_frag.groupby("name_bin").agg({"count": "sum"}).reset_index()
        )
        bedgraph_bins = self.binner.bins.merge(
            count_aggregated, left_on="name", right_on="name_bin"
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
                    df["count_overfrac_norm"] / self.n_cis_interactions
                )
                * scale_factor,
            )
        )

        count_aggregated = (
            bedgraph_frag.groupby("name_bin")
            .agg({"count_overfrac_n_interact_norm": "sum"})
            .reset_index()
        )

        bedgraph_bins = self.binner.bins.merge(
            count_aggregated, left_on="name", right_on="name_bin"
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

    def normalise_by_n_cis(self, reporter_distance=1e5, n_read_scale=1e6):

        if reporter_distance:
            n_reporters = (
                self.df.loc[lambda df: (df["chrom"] == self.capture_chrom)]
                .loc[lambda df: df["start"] >= (self.capture_start - reporter_distance)]
                .loc[lambda df: df["end"] <= (self.capture_end + reporter_distance)][
                    "score"
                ]
                .sum()
            )
        else:
            n_reporters = self.df.loc[lambda df: (df["chrom"] == self.capture_chrom)][
                "score"
            ].sum()

        return (self / n_read_scale) * n_reporters

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