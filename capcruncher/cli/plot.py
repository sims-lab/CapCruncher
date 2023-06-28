import os
import pathlib
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import capcruncher.api.plotting as cp
from loguru import logger

# def make_template(
#     files: list,
#     output_prefix: str,
#     design_matrix=None,
#     analysis_method="tiled",
#     viewpoint=None,
#     binsize=None,
#     genes=None,
# ):


#     tracks = []
#     tracks.append(cp.CCTrack(None, type="scale"))

#     if analysis_method != "tiled" and genes:
#         tracks.append(cp.CCTrack(genes, type="genes"))


#     ## Design matrix

#     # If a design matrix is provided, we need to group the files together
#     if design_matrix:
#         df_design = pd.read_csv(
#             design_matrix, sep=r",|\s+|\t", index_col="sample", engine="python"
#         )


#     # Deal with any files that need to be grouped together to form a summary dataframe
#     if design_matrix:

#         # Assuming design matrix has columns: sample  condition
#         df_design = pd.read_csv(
#             design_matrix, sep=r",|\s+|\t", index_col="sample", engine="python"
#         )
#         df_fnames = (
#             pd.Series(files).loc[lambda ser: ser.str.contains(".bigWig")].to_frame("fn")
#         )

#         # Extract sample name, normalisation and viewpoint from each file name
#         df_fnames["basename"] = df_fnames["fn"].apply(os.path.basename)
#         df_fnames = df_fnames.join(
#             df_fnames["basename"].str.extract(
#                 r"^(?P<samplename>.*?)\.(?P<normalisation>.*?)\.(?P<viewpoint>.*?)\.bigWig$"
#             )
#         )
#         df_fnames = df_fnames.set_index("samplename")
#         df_fnames = df_fnames.join(df_design["condition"])

#         # Need to ignore any subtraction bigWig files as these will interfere
#         df_fnames = df_fnames.loc[
#             lambda df: ~df["normalisation"].str.contains("-subtraction")
#         ]

#         # Get random colors for each plot
#         colors = [
#             matplotlib.colors.to_hex(c) for c in sns.palettes.hls_palette(n_colors=12)
#         ]

#         # Make bigwig collection and append to template
#         for ((condition, viewpoint), df), color in zip(
#             df_fnames.groupby(["condition", "viewpoint"]), colors
#         ):
#             fnames = df["fn"].to_list()
#             processed_files.update(fnames)
#             tracks[f"{condition}_{viewpoint}"] = bwc(file=fnames, color=color)

#     # Deal with the rest of the files
#     for fn in files:

#         if fn not in processed_files:

#             fn_base = os.path.basename(fn).replace(".gz", "")
#             fn_no_ext, ext = os.path.splitext(fn_base)
#             track_type = extensions_to_track_mapping.get(ext, None)

#             if track_type:
#                 tracks[fn_no_ext] = track_type(**{"file": fn})

#             else:
#                 raise ValueError(f"Track extension {ext} not supported at the moment")

#     tracks_for_output = {k: v._asdict() for k, v in tracks.items()}
#     with open(f"{output_prefix}.yml", "w") as w:
#         yaml.dump(tracks_for_output, w, sort_keys=False)


def plot(
    region: str,
    template: os.PathLike,
    output: str,
) -> None:
    """Plot a region using a template.

    Args:
        region (str): Genomic region to plot.
        template (os.PathLike): Path to template file.
        output (str): Path to output file.

    """

    fig = cp.CCFigure.from_toml(template)
    fig.save(region, output=output)
