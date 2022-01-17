import os
import sys
from collections import namedtuple, OrderedDict
import matplotlib
import pandas as pd
import yaml
import seaborn as sns

import coolbox.api as cb
from capcruncher.tools.plotting import (
    ScaleBar,
    SimpleBed,
    CCBigWig,
    CCBigWigCollection,
    CCMatrix,
    XAxisGenomic,
)


def make_template(
    files: list,
    output_prefix: str,
    design_matrix=None,
    analysis_method="tiled",
    viewpoint=None,
    binsize=None,
):

    bed = namedtuple(
        "Bed", field_names=["file", "color", "type"], defaults=["black", "bed"]
    )
    bw = namedtuple(
        "BigWig",
        field_names=["file", "color", "type", "style", "min_value", "max_value"],
        defaults=["blue", "bigWig", "fragment", "auto", "auto"],
    )
    bwc = namedtuple(
        "BigWigCollection",
        field_names=[
            "file",
            "color",
            "type",
            "smooth_window",
            "min_value",
            "max_value",
        ],
        defaults=["blue", "bigWigCollection", 2001, "auto", "auto"],
    )
    heatmap = namedtuple(
        "Heatmap",
        field_names=[
            "file",
            "viewpoint",
            "binsize",
            "assay",
            "normalization",
            "remove_viewpoint",
            "type",
            "color",
            "style",
            "min_value",
            "max_value",
            "transform",
            "orientation",
        ],
        defaults=[
            viewpoint if viewpoint else "VIEWPOINT",
            int(binsize) if binsize else "BINSIZE",
            analysis_method,
            "icen",
            False,
            "heatmap",
            "jet",
            "triangular",
            "auto",
            "auto",
            "no",
            "normal",
        ],
    )
    genes = namedtuple(
        "genes",
        field_names=["file", "type", "style", "gene_style", "color", "display", "fontsize"],
        defaults=["genes", "gene", "normal", "bed_rgb", "stacked", 7.5],
    )

    extensions_to_track_mapping = {
        "genes.bed": genes,
        ".bed": bed,
        ".bigWig": bw,
        ".hdf5": heatmap,
    }

    tracks = OrderedDict()
    processed_files = set()

    # Find any gene specifiers as these need to go first
    for fn in files:
        if "genes.bed" in fn:
            tracks[os.path.basename(fn)] = genes(file=fn)
            processed_files.add(fn)

    # Deal with any files that need to be grouped together to form a summary dataframe
    if design_matrix:

        # Assuming design matrix has columns: sample  condition
        df_design = pd.read_csv(design_matrix, sep=r",|\s+|\t", index_col="sample", engine="python")
        df_fnames = (
            pd.Series(files).loc[lambda ser: ser.str.contains(".bigWig")].to_frame("fn")
        )

        # Extract sample name, normalisation and viewpoint from each file name
        df_fnames["basename"] = df_fnames["fn"].apply(os.path.basename)
        df_fnames = df_fnames.join(
            df_fnames["basename"].str.extract(
                r"^(?P<samplename>.*?)\.(?P<normalisation>.*?)\.(?P<viewpoint>.*?)\.bigWig$"
            )
        )
        df_fnames = df_fnames.set_index("samplename")
        df_fnames = df_fnames.join(df_design["condition"])

        # Need to ignore any subtraction bigWig files as these will interfere
        df_fnames = df_fnames.loc[
            lambda df: ~df["normalisation"].str.contains("-subtraction")
        ]

        # Get random colors for each plot
        colors = [
            matplotlib.colors.to_hex(c) for c in sns.palettes.hls_palette(n_colors=12)
        ]

        # Make bigwig collection and append to template
        for ((condition, viewpoint), df), color in zip(
            df_fnames.groupby(["condition", "viewpoint"]), colors
        ):
            fnames = df["fn"].to_list()
            processed_files.update(fnames)
            tracks[f"{condition}_{viewpoint}"] = bwc(file=fnames, color=color)

    # Deal with the rest of the files
    for fn in files:

        if not fn in processed_files:

            fn_base = os.path.basename(fn).replace(".gz", "")
            fn_no_ext, ext = os.path.splitext(fn_base)
            track_type = extensions_to_track_mapping.get(ext, None)

            if track_type:
                tracks[fn_no_ext] = track_type(**{"file": fn})

            else:
                raise ValueError(f"Track extension {ext} not supported at the moment")

    tracks_for_output = {k: v._asdict() for k, v in tracks.items()}
    with open(f"{output_prefix}.yml", "w") as w:
        yaml.dump(tracks_for_output, w, sort_keys=False)


def plot_reporters(
    region: str, config: os.PathLike, output: str, x_axis=False, no_scale_bar=False,
):

    track_type_to_track_class_mapping = {
        "bigWig": CCBigWig,
        "bigWigCollection": CCBigWigCollection,
        "heatmap": CCMatrix,
        "bed": SimpleBed,
        "genes": cb.BED,
    }

    chrom = region.split(":")[0]
    start, end = [int(x) for x in region.split(":")[1].replace(",", "").split("-")]

    with open(config, "r") as r:
        tracks = yaml.safe_load(r)

    # Perform the plotting
    frame = cb.Frame()

    if not no_scale_bar:
        frame.add_track(ScaleBar())
    else:
        # Just add a spacer if no scale bar
        frame.add_track(cb.Spacer())

    for ii, (track_name, track_details) in enumerate(tracks.items()):

        track_type = track_details["type"]
        track_class = track_type_to_track_class_mapping[track_type]

        if track_class is cb.BED:
            # Don't give the bed file a title
            track = track_class(**track_details, title=" ")
        else:
            track = track_class(**track_details, title=track_name)

        frame.add_track(track)
        frame.add_track(cb.Spacer())

    if x_axis:
        frame.add_track(XAxisGenomic())

    figure = frame.plot(chrom, start, end)
    figure.savefig(output)
