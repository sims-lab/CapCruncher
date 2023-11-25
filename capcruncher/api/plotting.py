import functools
import math
import os
import pathlib
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

import coolbox.api as cb
import cooler.api as cooler
import iced
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
import tqdm
from coolbox.api import GenomeRange
from coolbox.core.track import Track
from coolbox.utilities import get_coverage_stack, get_feature_stack
from coolbox.utilities.genome import GenomeRange
from loguru import logger
from matplotlib import cm, colors, transforms
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from pybedtools import BedTool

import capcruncher.api as cc


class CCMatrix(cb.Cool):
    def __init__(
        self,
        file: os.PathLike,
        binsize: 5000,
        viewpoint: str,
        remove_viewpoint=False,
        **kwargs,
    ):
        self.binsize = binsize
        self.viewpoint = viewpoint
        self.remove_viewpoint = remove_viewpoint
        self.properties = dict()
        self.properties.update(kwargs)
        self.properties["name"] = f"CCMatrix.{self.properties.get('title')}"
        super(CCMatrix, self).__init__(file, **kwargs)
        # Need to override the coolbox default if we need a cmap to be set
        self.properties["color"] = kwargs.get("color", self.properties["color"])

        # Override the defaults
        self.properties["balance"] = "no"

        if not self._cooler_store_has_binsize:
            raise ValueError(
                f"Viewpoint {viewpoint} or resolution {binsize} not found in supplied file."
            )

        self.cooler = cooler.Cooler(f"{file}::{viewpoint}/resolutions/{binsize}")
        self.capture_bins = self.cooler.info["metadata"]["viewpoint_bins"]

    def _cooler_store_has_binsize(self):
        clrs = cooler.fileops.list_coolers(self.file)
        expected_path = f"{self.viewpoint}/resolutions/{self.binsize}"

        if expected_path in clrs:
            return True

    def get_matrix(self, coordinates, field="count"):
        matrix = self.cooler.matrix(field=field, balance=False).fetch(coordinates)

        offset = self.cooler.offset(coordinates)
        capture_bins = [(bin - offset) for bin in self.capture_bins]

        if self.remove_viewpoint:
            matrix[capture_bins, :] = 0
            matrix[:, capture_bins] = 0

        return matrix

    def get_matrix_normalised(
        self, coordinates, normalization_method=None, **normalisation_kwargs
    ):
        methods_stored = {
            "n_interactions": "count_n_interactions_norm",
            "n_rf_n_interactions": "count_n_rf_n_interactions_norm",
        }

        if normalization_method == "raw":
            matrix_normalised = self.get_matrix(coordinates)

        elif normalization_method in methods_stored:
            matrix_normalised = self.get_matrix(
                coordinates, field=methods_stored[normalization_method]
            )

        elif normalization_method == "ice":
            matrix = self.get_matrix(coordinates)
            matrix = np.nan_to_num(matrix)
            # matrix = iced.filter.filter_low_counts(matrix, percentage=0.04)
            matrix_normalised = iced.normalization.ICE_normalization(
                matrix, **normalisation_kwargs
            )  # Get iced matrix

        elif normalization_method == "icen_cis":
            matrix = self.get_matrix(coordinates)
            matrix = np.nan_to_num(matrix)
            matrix_ice = iced.normalization.ICE_normalization(
                matrix, **normalisation_kwargs
            )  # Get iced matrix
            matrix_normalised = (
                matrix_ice
                / int(self.cooler.info["metadata"]["n_cis_interactions"])
                * 1e6
            )  # Correct for number of interactions * 1e6

        elif normalization_method == "icen_scale":
            matrix = self.get_matrix(coordinates)
            matrix = np.nan_to_num(matrix)
            matrix_ice = iced.normalization.ICE_normalization(
                matrix, **normalisation_kwargs
            )  # Get iced matrix
            matrix_normalised = matrix_ice / self.properties["scaling_factor"]

        else:
            raise ValueError(
                f'Incorrect normalisation specified choose from: {" ".join(["raw", *methods_stored.keys(),"ice"])}'
            )

        return matrix_normalised

    def fetch_data(
        self, gr: cb.GenomeRange, gr2: cb.GenomeRange = None, **kwargs
    ) -> np.ndarray:
        norm = self.properties.get("normalization", "raw")
        matrix = self.get_matrix_normalised(
            f"{gr.chrom}:{gr.start}-{gr.end}", normalization_method=norm, **kwargs
        )
        return self.fill_zero_nan(matrix)

    def plot_matrix(self, gr: GenomeRange, gr2: GenomeRange = None):
        # Code taken and adapted from coolbox
        gr = GenomeRange(gr)

        if "JuiceBox" in self.properties["color"]:
            cmap = CCMatrix.get_juicebox_cmaps()[self.properties["color"]]
        else:
            cmap = cm.get_cmap(self.properties["color"])

        lowest = cmap(0)
        cmap.set_bad(lowest)
        cmap.set_under(lowest)

        ax = self.ax
        arr = self.matrix
        c_min, c_max = self.matrix_val_range

        if self.properties["max_value"] == "auto":
            matrix_triu = np.triu(self.matrix)
            c_max = np.percentile(matrix_triu, 98)

        if gr2 is None and self.style == self.STYLE_TRIANGULAR:
            # triangular style
            scale_r = 1 / math.sqrt(2)
            r_len = gr.end - gr.start
            # Rotate image using Affine2D, reference:
            #     https://stackoverflow.com/a/50920567/8500469

            tr = (
                transforms.Affine2D()
                .translate(-gr.start, -gr.start)
                .rotate_deg_around(0, 0, 45)
                .scale(scale_r)
                .translate(gr.start + r_len / 2, -r_len / 2)
            )

            img = ax.matshow(
                arr,
                cmap=cmap,
                transform=tr + ax.transData,
                extent=(gr.start, gr.end, gr.start, gr.end),
                aspect="auto",
                interpolation="none",
            )

        elif gr2 is None and self.style == self.STYLE_WINDOW:
            # window style
            # exist in HicMatBase
            fgr = self.fetched_gr
            scale_factor = fgr.length / gr.length
            scale_r = scale_factor / math.sqrt(2)
            length_dialog = gr.length * scale_factor
            delta_x = length_dialog * (gr.start - fgr.start) / fgr.length
            delta_x = length_dialog / 2 - delta_x
            tr = (
                transforms.Affine2D()
                .translate(-gr.start, -gr.start)
                .rotate_deg_around(0, 0, 45)
                .scale(scale_r)
                .translate(gr.start + delta_x, -fgr.length / 2)
            )
            img = ax.matshow(
                arr,
                cmap=cmap,
                transform=tr + ax.transData,
                extent=(gr.start, gr.end, gr.start, gr.end),
                aspect="auto",
            )
        else:
            if gr2 is None:
                gr2 = gr
            # matrix style
            img = ax.matshow(
                arr,
                cmap=cmap,
                extent=(gr.start, gr.end, gr2.end, gr2.start),
                aspect="auto",
            )

        if self.norm == "log":
            img.set_norm(colors.LogNorm(vmin=c_min, vmax=c_max))
        else:
            img.set_norm(colors.Normalize(vmin=c_min, vmax=c_max))

        return img

    @staticmethod
    def get_juicebox_cmaps():
        JuiceBoxLikeColor = LinearSegmentedColormap.from_list(
            "interaction", ["#FFFFFF", "#FFDFDF", "#FF7575", "#FF2626", "#F70000"]
        )
        JuiceBoxLikeColor.set_bad("white")
        JuiceBoxLikeColor.set_under("white")
        JuiceBoxLikeColor2 = LinearSegmentedColormap.from_list(
            "interaction", ["#FFFFFF", "#FFDFAF", "#FF7555", "#FF2600", "#F70000"]
        )
        JuiceBoxLikeColor2.set_bad("white")
        JuiceBoxLikeColor2.set_under("white")

        return {
            "JuiceBoxLike": JuiceBoxLikeColor,
            "JuiceBoxLike2": JuiceBoxLikeColor2,
        }


class CCBigWig(cb.BigWig):
    def __init__(self, file, **kwargs):
        self.file = file
        self.coverages = []

        super(CCBigWig, self).__init__(file, **kwargs)

    def fetch_data(self, gr, **kwargs):
        if not self.properties["style"] == "fragment":
            data = super(CCBigWig, self).fetch_data(gr, **kwargs)
        else:
            data = self.bw.fetch_intervals(gr.chrom, gr.start, gr.end)

        return data

    def plot_fragments(self, ax, gr, **kwargs):
        data = self.fetch_data(gr, **kwargs)
        _alpha = self.properties.get("alpha", 1.0)
        _threshold = self.properties.get("threshold", 0)
        _offset = gr.start
        bp_proportion = 1 / (data["end"].max() - data["start"].min())

        for row in data.itertuples():
            pg = Polygon(
                [
                    (row.start, 0),
                    (row.start, row.value),
                    (row.end, row.value),
                    (row.end, 0),
                ],
                color=self.properties["color"],
            )
            ax.add_patch(pg)

        ax.set_ylim(0, data["value"].max())
        ymin, ymax = self.adjust_plot(ax, gr)
        self.plot_data_range(ax, ymin, ymax, self.properties["data_range_style"], gr)
        self.plot_label()

    def plot(self, ax, gr, **kwargs):
        if not self.properties["style"] == "fragment":
            super(CCBigWig, self).plot(ax, gr, **kwargs)
        else:
            self.plot_fragments(ax, gr, **kwargs)


class CCBigWigCollection(Track):
    DEFAULT_PROPERTIES = {
        "style": "line",
        "fmt": "-",
        "line_width": 2.0,
        "size": 10,
        "color": "#a6cee3",
        "threshold_color": "#ff9c9c",
        "threshold": "inf",
        "cmap": "bwr",
        "orientation": None,
        "data_range_style": "y-axis",
        "min_value": "auto",
        "max_value": "auto",
    }

    def __init__(self, file: list, exclusions: str = None, **kwargs):
        self.file_names = file
        self.exclusions = exclusions
        self.bws = [cb.BigWig(str(fn)) for fn in file]
        self.properties = {"files": self.file_names}
        self.properties.update(CCBigWigCollection.DEFAULT_PROPERTIES.copy())
        self.properties.update(kwargs)
        self.properties["name"] = f"BigWigCollection.{self.properties.get('title')}"
        super(CCBigWigCollection, self).__init__(**self.properties)

        self.coverages = []

        # load features from global feature stack
        features_stack = get_feature_stack()
        for features in features_stack:
            self.properties.update(features.properties)

        # load coverages from global coverages stack
        coverage_stack = get_coverage_stack()
        for coverage in coverage_stack:
            self.coverages.append(coverage)

    def fetch_data(self, gr, **kwargs):
        datasets = [
            bw.bw.fetch_intervals(gr.chrom, gr.start, gr.end)
            .set_index(["chrom", "start", "end"])
            .rename(columns={"value": os.path.basename(bw.properties["file"])})
            for bw in self.bws
        ]
        df = datasets[0].join(datasets[1:])
        df_summary = df.assign(mean=df.mean(axis=1), sem=df.sem(axis=1)).reset_index()

        intervals_to_bp = []
        for interval in df_summary.itertuples():
            interval_len = interval.end - interval.start

            interval_positions = np.arange(interval_len) + interval.start
            scores_mean = np.repeat(interval.mean, interval_len)
            scores_sem = np.repeat(interval.sem, interval_len)

            intervals_to_bp.append(
                np.vstack([interval_positions, scores_mean, scores_sem]).T
            )

        df_intervals = pd.DataFrame(
            np.concatenate(intervals_to_bp), columns=["bp", "mean", "sem"]
        )

        if self.exclusions:
            df_intervals = pd.concat(
                [df_intervals, self.fetch_exluded_regions(gr)]
            ).sort_values("bp")

        if self.properties.get("smooth_window"):
            from scipy.signal import savgol_filter

            df_intervals["mean_smoothed"] = savgol_filter(
                df_intervals["mean"],
                window_length=self.properties.get("smooth_window", 1001),
                polyorder=self.properties.get("polyorder", 1),
            )

        return df_intervals

    def fetch_exluded_regions(self, gr):
        excluded_tabix = BedTool(self.exclusions).tabix(force=True)
        df_excluded = excluded_tabix.tabix_intervals(
            f"{gr.chrom}:{gr.start}-{gr.end}"
        ).to_dataframe()

        intervals_to_bp = []
        for interval in df_excluded.itertuples():
            interval_len = interval.end - interval.start

            interval_positions = np.arange(interval_len) + interval.start
            scores_nan = np.repeat(np.nan, interval_len)
            intervals_to_bp.append(interval_positions)

        df_intervals = pd.Series(np.concatenate(intervals_to_bp)).to_frame("bp")
        df_intervals["mean"] = np.nan
        df_intervals["sem"] = np.nan

        return df_intervals

    def plot(self, ax, gr, **kwargs):
        data = self.fetch_data(gr, **kwargs)

        line_width = self.properties.get("line_width", 1)
        color = self.properties.get("color", "blue")
        alpha = self.properties.get("alpha", 0.2)
        downsample = self.properties.get("downsample", 0)

        if downsample:
            rows = np.arange(0, data.shape[0], downsample)
            data = data.iloc[rows]

        if self.properties.get("smooth_window"):
            scores = data["mean_smoothed"]
        else:
            scores = data["mean"]

        ax.fill_between(
            data["bp"],
            scores - data["sem"],
            scores + data["sem"],
            alpha=alpha,
            color=color,
            zorder=0,
        )

        ax.plot(
            data["bp"],
            scores,
            color=color,
            zorder=1,
        )

        min_val = self.properties.get("min_value")
        max_val = self.properties.get("max_value")

        ymin = round(scores.min()) if min_val == "auto" else min_val
        ymax = round(scores.max() + data["sem"].max()) if max_val == "auto" else max_val

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(ymin, ymax)

        self.plot_data_range(ax, ymin, ymax, self.properties["data_range_style"], gr)
        self.plot_label()

    def plot_data_range(self, ax, ymin, ymax, data_range_style, gr: cb.GenomeRange):
        if data_range_style == "text":
            self.plot_text_range(ax, ymin, ymax, gr)
        else:  # 'y-axis' style
            try:
                y_ax = self.y_ax
                self.plot_yaxis_range(ax, y_ax)
            except AttributeError:
                self.plot_data_range(ax, ymin, ymax, "text", gr)

    def plot_yaxis_range(self, plot_axis, y_ax):
        # """
        # Plot the scale of the y axis with respect to the plot_axis
        # plot something that looks like this:
        # ymax ┐
        #      │
        #      │
        # ymin ┘
        # Parameters
        # ----------
        # plot_axis : matplotlib.axes.Axes
        #     Main plot axis.
        # y_ax : matplotlib.axes.Axes
        #     Axis to use to plot the scale
        # """

        if (
            "show_data_range" in self.properties
            and self.properties["show_data_range"] == "no"
        ):
            return

        def value_to_str(value):
            if value % 1 == 0:
                return str(int(value))
            else:
                return f"{value:.4f}" if value < 0.01 else f"{value:.2f}"

        ymin, ymax = plot_axis.get_ylim()

        ymax_str = value_to_str(ymax)
        ymin_str = value_to_str(ymin)
        x_pos = [0, 0.5, 0.5, 0]
        y_pos = [0.01, 0.01, 0.99, 0.99]
        y_ax.plot(x_pos, y_pos, color="black", linewidth=1, transform=y_ax.transAxes)
        y_ax.text(
            -0.2,
            -0.01,
            ymin_str,
            verticalalignment="bottom",
            horizontalalignment="right",
            transform=y_ax.transAxes,
        )
        y_ax.text(
            -0.2,
            1,
            ymax_str,
            verticalalignment="top",
            horizontalalignment="right",
            transform=y_ax.transAxes,
        )
        y_ax.patch.set_visible(False)

    def plot_text_range(self, ax, ymin, ymax, gr: cb.GenomeRange):
        ydelta = ymax - ymin

        # set min max
        def format_lim(lim):
            return int(lim) if float(lim) % 1 == 0 else f"{lim:.2f}"

        ymax_print = format_lim(ymax)
        ymin_print = format_lim(ymin)
        small_x = 0.01 * gr.length
        # by default show the data range
        ax.text(
            gr.start - small_x,
            ymax - ydelta * 0.2,
            f"[ {ymin_print} ~ {ymax_print} ]",
            horizontalalignment="left",
            verticalalignment="top",
        )


class ScaleBar(Track):
    def __init__(self, **kwargs):
        self.properties = dict()
        self.properties["name"] = "Scale"
        self.properties.update(kwargs)
        super(ScaleBar, self).__init__()

    def fetch_data(self, **kwargs):
        pass

    def get_appropriate_scale(self, length):
        if length <= 1e3:
            scale = 1e2
        elif 1e3 < length < 1e4:
            scale = 1e3
        elif 1e4 < length < 1e5:
            scale = 1e4
        elif 1e5 < length < 1e6:
            scale = 1e5
        elif 1e6 < length < 1e7:
            scale = 1e6
        elif 1e7 < length < 1e8:
            scale = 1e8

        return scale

    def plot(self, ax, gr, **kwargs):
        position = self.properties.get("position", "left")
        y_midpoint = 0.5

        if self.properties.get("scale_distance"):
            scale_distance = self.properties["scale_distance"]

        else:
            scale_distance = self.get_appropriate_scale(gr.end - gr.start)

        # Determine x start and end based on position
        if position == "left":
            x0 = gr.start
            x1 = x0 + scale_distance
        elif position == "right":
            x0 = gr.end - scale_distance
            x1 = gr.end
        else:
            raise ValueError('Position can only be "left" or "right"')

        # Plot scale bar
        ax.plot([x0, x1], [y_midpoint, y_midpoint], color="black")
        ax.plot([x0, x0], [y_midpoint - 0.1, y_midpoint + 0.1], color="black", lw=1)
        ax.plot([x1, x1], [y_midpoint - 0.1, y_midpoint + 0.1], color="black", lw=1)

        # Add annotation
        from capcruncher.utils import get_human_readable_number_of_bp

        scale_distance_human_readable = get_human_readable_number_of_bp(scale_distance)

        ax.text(
            (x0 + (scale_distance / 2)),
            y_midpoint - 0.2,
            scale_distance_human_readable,
            ha="center",
            va="center",
        )

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)


class CCSimpleBed(cb.BED):
    def __init__(self, file: str, **kwargs):
        self.file = file
        self.properties = dict()
        self.properties["name"] = "BlockBed"
        self.properties.update(kwargs)

    def fetch_data(self, gr):
        bt = BedTool(self.file)
        bt_tabix = bt.tabix(force=True)

        return bt_tabix.tabix_intervals(
            f"{gr.chrom}:{gr.start}-{gr.end}"
        ).to_dataframe()

    def plot(self, ax, gr, **kwargs):
        data = self.fetch_data(gr)
        y_midpoint = 0.5

        for interval in data.itertuples():
            pg = Polygon(
                [
                    (interval.start, y_midpoint - 0.1),
                    (interval.start, y_midpoint + 0.1),
                    (interval.end, y_midpoint + 0.1),
                    (interval.end, y_midpoint - 0.1),
                ],
                color=self.properties.get("color", "black"),
            )

            if hasattr(interval, "name") and not self.properties.get("no_annotation"):
                interval_midpoint = interval.start + (
                    (interval.end - interval.start) / 2
                )
                ax.text(
                    interval_midpoint,
                    y_midpoint - 0.1,
                    interval.name,
                    ha="center",
                    va="center",
                )

            ax.add_patch(pg)

        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)


class CCXAxisGenomic(cb.XAxis):
    def __init__(self, **kwargs):
        super(CCXAxisGenomic, self).__init__()
        self.properties.update(kwargs)

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax

        ax.set_xlim(gr.start, gr.end)
        ticks = np.linspace(gr.start, gr.end, 10)
        labels = [f"{x:,.0f}" for x in ticks]

        ax.axis["x"] = ax.new_floating_axis(0, 0.5)
        ax.axis["x"].axis.set_ticks(ticks)
        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis["x"].axis.set_tick_params(which="minor", bottom="on")

        ax.axis["x"].major_ticklabels.set(size=10)

        if "where" in self.properties and self.properties["where"] == "top":
            ax.axis["x"].set_axis_direction("top")


# class CCMatrix(cb.Cool):
#     def __init__(
#         self,
#         file: os.PathLike,
#         binsize: 5000,
#         viewpoint: str,
#         remove_viewpoint=False,
#         **kwargs,
#     ):
#         self.binsize = binsize
#         self.viewpoint = viewpoint
#         self.remove_viewpoint = remove_viewpoint
#         self.properties = dict()
#         self.properties.update(kwargs)
#         self.properties["name"] = f"CCMatrix.{self.properties.get('title')}"
#         super(CCMatrix, self).__init__(file, **kwargs)
#         # Need to override the coolbox default if we need a cmap to be set
#         self.properties["color"] = kwargs.get("color", self.properties["color"])

#         # Override the defaults
#         self.properties["balance"] = "no"

#         if not self._cooler_store_has_binsize:
#             raise ValueError(
#                 f"Viewpoint {viewpoint} or resolution {binsize} not found in supplied file."
#             )

#         self.cooler = cooler.Cooler(f"{file}::{viewpoint}/resolutions/{binsize}")
#         self.capture_bins = self.cooler.info["metadata"]["viewpoint_bins"]

#     def _cooler_store_has_binsize(self):
#         clrs = cooler.fileops.list_coolers(self.file)
#         expected_path = f"{self.viewpoint}/resolutions/{self.binsize}"

#         if expected_path in clrs:
#             return True

#     def get_matrix(self, coordinates, field="count"):
#         matrix = self.cooler.matrix(field=field, balance=False).fetch(coordinates)

#         offset = self.cooler.offset(coordinates)
#         capture_bins = [(bin - offset) for bin in self.capture_bins]

#         if self.remove_viewpoint:
#             matrix[capture_bins, :] = 0
#             matrix[:, capture_bins] = 0

#         return matrix

#     def get_matrix_normalised(
#         self, coordinates, normalization_method=None, **normalisation_kwargs
#     ):
#         methods_stored = {
#             "n_interactions": "count_n_interactions_norm",
#             "n_rf_n_interactions": "count_n_rf_n_interactions_norm",
#         }

#         if normalization_method == "raw":
#             matrix_normalised = self.get_matrix(coordinates)

#         elif normalization_method in methods_stored:
#             matrix_normalised = self.get_matrix(
#                 coordinates, field=methods_stored[normalization_method]
#             )

#         elif normalization_method == "ice":
#             matrix = self.get_matrix(coordinates)
#             matrix = np.nan_to_num(matrix)
#             # matrix = iced.filter.filter_low_counts(matrix, percentage=0.04)
#             matrix_normalised = iced.normalization.ICE_normalization(
#                 matrix, **normalisation_kwargs
#             )  # Get iced matrix

#         elif normalization_method == "icen_cis":
#             matrix = self.get_matrix(coordinates)
#             matrix = np.nan_to_num(matrix)
#             matrix_ice = iced.normalization.ICE_normalization(
#                 matrix, **normalisation_kwargs
#             )  # Get iced matrix
#             matrix_normalised = (
#                 matrix_ice
#                 / int(self.cooler.info["metadata"]["n_cis_interactions"])
#                 * 1e6
#             )  # Correct for number of interactions * 1e6

#         elif normalization_method == "icen_scale":
#             matrix = self.get_matrix(coordinates)
#             matrix = np.nan_to_num(matrix)
#             matrix_ice = iced.normalization.ICE_normalization(
#                 matrix, **normalisation_kwargs
#             )  # Get iced matrix
#             matrix_normalised = matrix_ice / self.properties["scaling_factor"]

#         else:
#             raise ValueError(
#                 f'Incorrect normalisation specified choose from: {" ".join(["raw", *methods_stored.keys(),"ice"])}'
#             )

#         return matrix_normalised

#     def fetch_data(
#         self, gr: cb.GenomeRange, gr2: cb.GenomeRange = None, **kwargs
#     ) -> np.ndarray:
#         norm = self.properties.get("normalization", "raw")
#         matrix = self.get_matrix_normalised(
#             f"{gr.chrom}:{gr.start}-{gr.end}", normalization_method=norm, **kwargs
#         )
#         return self.fill_zero_nan(matrix)

#     def plot_matrix(self, gr: GenomeRange, gr2: GenomeRange = None):
#         # Code taken and adapted from coolbox
#         gr = GenomeRange(gr)

#         if "JuiceBox" in self.properties["color"]:
#             cmap = CCMatrix.get_juicebox_cmaps()[self.properties["color"]]
#         else:
#             cmap = cm.get_cmap(self.properties["color"])

#         lowest = cmap(0)
#         cmap.set_bad(lowest)
#         cmap.set_under(lowest)

#         ax = self.ax
#         arr = self.matrix
#         c_min, c_max = self.matrix_val_range

#         if self.properties["max_value"] == "auto":
#             matrix_triu = np.triu(self.matrix)
#             c_max = np.percentile(matrix_triu, 98)

#         if gr2 is None and self.style == self.STYLE_TRIANGULAR:
#             # triangular style
#             scale_r = 1 / math.sqrt(2)
#             r_len = gr.end - gr.start
#             # Rotate image using Affine2D, reference:
#             #     https://stackoverflow.com/a/50920567/8500469

#             tr = (
#                 transforms.Affine2D()
#                 .translate(-gr.start, -gr.start)
#                 .rotate_deg_around(0, 0, 45)
#                 .scale(scale_r)
#                 .translate(gr.start + r_len / 2, -r_len / 2)
#             )

#             img = ax.matshow(
#                 arr,
#                 cmap=cmap,
#                 transform=tr + ax.transData,
#                 extent=(gr.start, gr.end, gr.start, gr.end),
#                 aspect="auto",
#                 interpolation="none",
#             )

#         elif gr2 is None and self.style == self.STYLE_WINDOW:
#             # window style
#             # exist in HicMatBase
#             fgr = self.fetched_gr
#             scale_factor = fgr.length / gr.length
#             scale_r = scale_factor / math.sqrt(2)
#             length_dialog = gr.length * scale_factor
#             delta_x = length_dialog * (gr.start - fgr.start) / fgr.length
#             delta_x = length_dialog / 2 - delta_x
#             tr = (
#                 transforms.Affine2D()
#                 .translate(-gr.start, -gr.start)
#                 .rotate_deg_around(0, 0, 45)
#                 .scale(scale_r)
#                 .translate(gr.start + delta_x, -fgr.length / 2)
#             )
#             img = ax.matshow(
#                 arr,
#                 cmap=cmap,
#                 transform=tr + ax.transData,
#                 extent=(gr.start, gr.end, gr.start, gr.end),
#                 aspect="auto",
#             )
#         else:
#             if gr2 is None:
#                 gr2 = gr
#             # matrix style
#             img = ax.matshow(
#                 arr,
#                 cmap=cmap,
#                 extent=(gr.start, gr.end, gr2.end, gr2.start),
#                 aspect="auto",
#             )

#         if self.norm == "log":
#             img.set_norm(colors.LogNorm(vmin=c_min, vmax=c_max))
#         else:
#             img.set_norm(colors.Normalize(vmin=c_min, vmax=c_max))

#         return img

#     @staticmethod
#     def get_juicebox_cmaps():
#         JuiceBoxLikeColor = LinearSegmentedColormap.from_list(
#             "interaction", ["#FFFFFF", "#FFDFDF", "#FF7575", "#FF2626", "#F70000"]
#         )
#         JuiceBoxLikeColor.set_bad("white")
#         JuiceBoxLikeColor.set_under("white")
#         JuiceBoxLikeColor2 = LinearSegmentedColormap.from_list(
#             "interaction", ["#FFFFFF", "#FFDFAF", "#FF7555", "#FF2600", "#F70000"]
#         )
#         JuiceBoxLikeColor2.set_bad("white")
#         JuiceBoxLikeColor2.set_under("white")

#         return {
#             "JuiceBoxLike": JuiceBoxLikeColor,
#             "JuiceBoxLike2": JuiceBoxLikeColor2,
#         }


class CCAverageMatrix(CCMatrix):
    def __init__(
        self,
        matricies: List[CCMatrix],
        aggregation: Literal["sum", "mean", "median"] = "mean",
        **kwargs,
    ):
        self.matricies = matricies
        self.aggregation = aggregation
        self.properties = matricies[0].properties
        self.properties.update(kwargs)
        self.properties["name"] = f"CCMatrix.{self.properties.get('title')}"

        # Need to override the coolbox default if we need a cmap to be set
        self.properties["color"] = kwargs.get("color", self.properties["color"])

        # Override the defaults
        self.properties["balance"] = "no"

    @functools.cache
    def fetch_data(self, gr: cb.GenomeRange, gr2=None, **kwargs):
        data = []
        for matrix in tqdm.tqdm(self.matricies):
            data.append(matrix.fetch_data(gr, **kwargs))

        try:
            func_agg = getattr(np, self.aggregation)
        except AttributeError:
            raise ValueError(
                f"Aggregation function {self.aggregation} not found in numpy"
            )

        # Aggregate the list of matricies into a single matrix
        data = func_agg(data, axis=0)

        self.fetched_gr = gr
        self.fetched_gr2 = gr2
        return data

    def plot(self, ax, gr: GenomeRange, **kwargs):
        self.ax = ax
        # fetch processed plot_data
        self.matrix = self.fetch_data(gr, **kwargs)
        # plot matrix
        img = self.plot_matrix(gr, kwargs.get("gr2"))
        self.adjust_figure(gr, kwargs.get("gr2"))
        self.draw_colorbar(img)
        self.plot_label()


class CCTrack:
    """
    Provides a wrapper around tracks to provide a consistent interface

    Args:
        file (os.PathLike): Path to file to plot
        file_type (Literal["heatmap", "bigwig", "bigwig_summary", "scale", "bed", "xaxis", "genes", "spacer"], optional): Type of file to plot. Defaults to None.
        **kwargs: Additional arguments to pass to the track
    """

    def __init__(
        self,
        file,
        file_type: Literal[
            "heatmap",
            "heatmap_summary",
            "bigwig",
            "bigwig_summary",
            "scale",
            "bed",
            "xaxis",
            "genes",
            "spacer",
        ] = None,
        **kwargs,
    ):
        self.file = file
        self.properties = dict()
        self.properties.update(kwargs)

        if file_type:
            self.properties["type"] = file_type
        elif self.properties.get("type"):
            pass
        else:
            raise ValueError(
                "Please specify file_type as one of: heatmap, bigwig, bigwig_summary, scale, bed, xaxis, genes, spacer"
            )

    def get_track(self):
        match self.properties.get("type"):  # noqa
            case "heatmap":
                assert (
                    "binsize" in self.properties
                ), "Binsize must be specified for heatmap track (e.g. binsize=5000)"
                return CCMatrix(self.file, **self.properties)
            case "heatmap_summary":
                assert (
                    "binsize" in self.properties
                ), "Binsize must be specified for heatmap track (e.g. binsize=5000)"

                matricies = []
                for matrix in self.file:
                    matricies.append(CCMatrix(matrix, **self.properties))

                return CCAverageMatrix(matricies, **self.properties)
            case "bigwig":
                if self.properties.get("overlay"):
                    return cb.BigWigCoverage(self.file, **self.properties)
                else:
                    return CCBigWig(self.file, **self.properties)
            case "bigwig_summary":
                return CCBigWigCollection(self.file, **self.properties)
            case "scale":
                return ScaleBar(**self.properties)
            case "bed":
                return CCSimpleBed(self.file, **self.properties)
            case "xaxis":
                return CCXAxisGenomic(**self.properties)
            case "genes":
                if self.properties.get("title"):
                    del self.properties["title"]
                return cb.BED(self.file, **self.properties)
            case "spacer":
                return cb.Spacer(**self.properties)
            case _:
                if getattr(cb, self.properties.get("type")):
                    return getattr(cb, self.properties.get("type"))(
                        self.file, **self.properties
                    )

                else:
                    raise ValueError(
                        f"Unknown track type {self.properties.get('type')}, select from: heatmap, bigwig, bigwig_summary, scale, bed, xaxis, genes, spacer"
                    )

    @property
    def path(self) -> str:
        if isinstance(self.file, (list, tuple, np.ndarray)):
            return [str(pathlib.Path(f).resolve()) for f in self.file]
        else:
            return str(pathlib.Path(self.file).resolve())

    def __repr__(self) -> str:
        return f"CCTrack({self.properties.get('title')}, {self.properties.get('type')})"


class CCFigure:
    """
    Generates a figure from a list of tracks

    Args:
        tracks (List[CCTrack], optional): List of tracks to plot. Defaults to None.
        auto_spacing (bool, optional): Automatically add a spacer track between each track. Defaults to False.
        **kwargs: Additional arguments to pass to the figure
    """

    def __init__(
        self, tracks: List[CCTrack] = None, auto_spacing: bool = False, **kwargs
    ) -> None:
        self.frame = cb.Frame()
        self.auto_spacing = auto_spacing
        self.properties = dict()
        self.properties.update(kwargs)

        if tracks:
            self.tracks = set(tracks)
            self.add_tracks(tracks)
        else:
            self.tracks = set()

    def add_track(self, track: CCTrack) -> None:
        """
        Add a track to the figure

        Args:
            track (CCTrack): Track to add
        """
        self.tracks.add(track)
        self.frame.add_track(track.get_track())

    def add_tracks(self, tracks: List[CCTrack]) -> None:
        """
        Add a list of tracks to the figure

        Args:
            tracks (List[CCTrack]): List of tracks to add
        """
        for track in tracks:
            if self.auto_spacing:
                spacer = CCTrack(None, file_type="spacer")
                self.add_track(spacer.get_track())

            self.add_track(track)

    def plot(
        self,
        gr: Union[str, GenomeRange],
        gr2: Union[str, GenomeRange] = None,
        show: bool = True,
        **kwargs,
    ) -> None:
        """
        Plot the figure

        Args:
            gr (Union[str, GenomeRange]): GenomeRange to plot
            gr2 (Union[str, GenomeRange], optional): Second GenomeRange to plot. Defaults to None.
            show (bool, optional): Show the figure. Defaults to True.
            **kwargs: Additional arguments to pass to the plot
        """

        if gr2:
            fig = self.frame.plot(gr, gr2, **kwargs)
        else:
            fig = self.frame.plot(gr, **kwargs)
        if show:
            fig.show()

        return fig

    def save(
        self,
        gr: Union[str, GenomeRange],
        gr2: Union[str, GenomeRange] = None,
        output: str = None,
        **kwargs,
    ) -> None:
        """
        Plots the figure and saves it to a file

        Args:
            gr (Union[str, GenomeRange]): GenomeRange to plot
            gr2 (Union[str, GenomeRange], optional): Second GenomeRange to plot. Defaults to None.
            output (str, optional): Path to save the figure to. Defaults to None.
            **kwargs: Additional arguments to pass to the plot
        """

        fig = self.plot(gr, gr2, show=False, **kwargs)
        if output:
            fig.savefig(output, dpi=300)
        else:
            fig.savefig(f"{gr.chrom}_{gr.start}_{gr.end}.png", dpi=300)

    @classmethod
    def from_toml(cls, toml_file: os.PathLike, **kwargs) -> "CCFigure":
        """
        Instantiate a CCFigure from a toml file

        Args:
            toml_file (os.PathLike): Path to toml file
            **kwargs: Additional arguments to pass to the figure
        """
        import toml

        with open(toml_file) as f:
            config = toml.load(f)

        tracks = []
        for track_name, attr in config.items():
            file = attr.pop("file") if attr.get("file") else None
            track_name = attr.pop("title") if attr.get("title") else track_name
            tracks.append(CCTrack(file, title=track_name, **attr))
        return cls(tracks, **kwargs)

    @classmethod
    def from_frame(cls, frame: cb.Frame, **kwargs) -> "CCFigure":
        """
        Instantiate a CCFigure from a coolbox Frame

        Args:
            frame (cb.Frame): coolbox Frame to instantiate from
            **kwargs: Additional arguments to pass to the figure
        """
        tracks = []
        for track in frame.tracks:
            tracks.append(CCTrack(track.properties["file"], **track.properties))

        return cls(tracks, **kwargs)

    def to_toml(self, output: str = None) -> Union[None, Dict[str, Any]]:
        """
        Save the CCFigure to a toml file

        Args:
            output (str, optional): Path to save the toml file to. Defaults to None.

        Returns:
            Union[None, Dict[str, Any]]: If output is not specified, returns a dict of the toml file

        """

        import toml
        from collections import OrderedDict

        def _get_n_tracks_of_type(config: Dict[str, Dict], track_type: str):
            return sum(1 for t in config.keys() if track_type in t)

        config = OrderedDict()
        for track in self.tracks:
            # Perform conversions for file-less tracks
            if track.properties.get("type") in ["spacer", "scale", "xaxis"]:
                track_type = track.properties.get("type")
                n = _get_n_tracks_of_type(config, track_type)
                config[f"{track_type} {n}"] = track.properties
                config[f"{track_type} {n}"]["file"] = None
            elif track.properties.get("type") == "genes":
                track_type = track.properties.get("type")
                n = _get_n_tracks_of_type(config, track_type)
                config[f"{track_type} {n}"] = track.properties
                config[f"{track_type} {n}"]["file"] = track.path
            else:
                config[track.properties["title"]] = track.properties
                config[track.properties["title"]]["file"] = track.path

        outfile = output if output else "config.toml"

        with open(outfile, "w") as f:
            config_str = toml.dumps(config)
            f.write(config_str)

        return config
