import os
import sys
from collections import namedtuple
import pandas as pd
import yaml
import glob

import click
import coolbox.api as cb
from capcruncher.tools.plotting import (
    ScaleBar,
    SimpleBed,
    CCBigWig,
    CCBigWigCollection,
    CCMatrix,
)
from capcruncher.cli.cli_reporters import cli


@cli.command()
@click.argument("files", nargs=-1)
@click.option("-o", "--output_prefix", default="template")
@click.option("-d", "--design_matrix")
def generate_plotting_template(
    files: list, output_prefix: str, design_matrix=None, analysis_method="tiled"
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
        field_names=["file", "color", "type", "smooth_window"],
        defaults=["blue", "bigWigCollection", 2001],
    )
    heatmap = namedtuple(
        "Heatmap",
        field_names=[
            "file",
            "viewpoint",
            "assay",
            "normalisation",
            "remove_viewpoint",
            "type",
            "color",
        ],
        defaults=[
            "VIEWPOINT",
            analysis_method,
            "ice",
            "remove_viewpoint",
            "heatmap",
            "red",
        ],
    )
    genes = namedtuple(
        "genes",
        field_names=["file", "type", "style", "gene_style", "color", "display"],
        defaults=["gene", "gene", "normal", "bed_rgb", "stacked"],
    )

    extensions_to_track_mapping = {
        ".bed": bed,
        ".bigWig": bw,
        ".hdf5": heatmap,
        "genes.bed": genes,
    }

    tracks = dict()
    processed_files = set()

    # Deal with any files that need to be grouped together to form a summary dataframe
    if design_matrix:
        # Assuming design matrix has columns: sample  condition
        df_design = pd.read_csv(design_matrix, sep="\t", index_col="sample")
        df_fnames = (
            pd.Series(files).loc[lambda ser: ser.str.contains(".bigWig")].to_frame("fn")
        )
        df_fnames["basename"] = df_fnames["fn"].apply(os.path.basename)
        df_fnames = df_fnames.join(
            df_fnames["basename"].str.extract(
                r"^(?P<samplename>.*?)\.(?P<normalisation>.*?)\.(?P<viewpoint>.*?)\.bigWig$"
            )
        )
        df_fnames = df_fnames.set_index("samplename")
        df_fnames = df_fnames.join(df_design["condition"])

        for (condition, viewpoint), df in df_fnames.groupby(["condition", "viewpoint"]):
            fnames = df["fn"].to_list()
            processed_files.update(fnames)
            tracks[f"{condition}_{viewpoint}"] = bwc(file=fnames)

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
        yaml.dump(tracks_for_output, w)


@cli.command()
def plot2():
    # def plot2(config: os.PathLike, output_prefix: str):

    region = "chr16:53951-278829"

    track_type_to_track_class_mapping = {
        "bigWig": CCBigWig,
        "bigWigCollection": CCBigWigCollection,
        "heatmap": CCMatrix,
        "bed": SimpleBed,
        "genes": cb.BED,
    }

    template = "template.yml"
    with open(template, "r") as r:
        tracks = yaml.load(r)

    # Perform the plotting
    frame = cb.Frame()
    frame.add_track(ScaleBar())

    for ii, (track_name, track_details) in enumerate(tracks.items()):

        track_type = track_details["type"]
        track_class = track_type_to_track_class_mapping[track_type]

        track = track_class(**track_details, title=track_name)

        frame.add_track(track)
        frame.add_track(cb.Spacer())

    figure = frame.plot(region)
    figure.savefig("plotting.png")
