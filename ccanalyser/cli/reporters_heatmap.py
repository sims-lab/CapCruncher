import os
import sys
from typing import Tuple, Union

import cooler
import pandas as pd
import numpy as np
import click

from ccanalyser.cli.cli_reporters import cli
from ccanalyser.utils import convert_interval_to_coords, format_coordinates
from ccanalyser.tools.plotting import CCMatrix


def plot_matrix(matrix, figsize=(10, 10), axis_labels=None, cmap=None, vmin=0, vmax=0):

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not vmax > 0:
        vmax = np.percentile(matrix, 95)
    
    if not vmin > 0:
        vmin = np.percentile(matrix, 5)

    fig, ax = plt.subplots()
    ax.imshow(matrix, vmax=vmax, vmin=vmin, cmap=cmap)
    ax.axis("off")

    return fig


@cli.command()
@click.argument("cooler_fn", required=True)
@click.option(
    "-c",
    "--coordinates",
    help="Coordinates in the format chr1:1000-2000 or a path to a .bed file with coordinates",
)
@click.option(
    "-r",
    "--resolution",
    help="Resolution at which to plot. Must be present within the cooler file.",
    required=True,
    multiple=True,
)
@click.option(
    "-n",
    "--capture_names",
    multiple=True,
    help="Capture names to plot. If not supplied will plot all",
)
@click.option(
    "--normalisation",
    help="Normalisation method for heatmap",
    default="n_interactions",
    type=click.Choice(["n_interactions", "n_rf_n_interactions", "ice"]),
)
@click.option("--cmap", help="Colour map to use for plotting", default="jet")
@click.option("--vmax", help="Vmax for plotting", default=0, type=click.FLOAT)
@click.option("--vmin", help="Vmin for plotting", default=0, type=click.FLOAT)
@click.option("-o", "--output_prefix", help="Output prefix for plot", default="")
@click.option("--remove_capture", help='Removes the capture probe bins from the matrix', is_flag=True, default=False)
def plot(
    cooler_fn: os.PathLike,
    coordinates: Union[str, os.PathLike],
    resolution: int,
    capture_names: Tuple = None,
    normalisation: str = None,
    cmap: str = "jet",
    vmax: float = 1,
    vmin: float = 0, 
    output_prefix: os.PathLike = "",
    remove_capture: bool = False,
):
    """
    Plots a heatmap of reporter interactions.

    Parses a HDF5 file containg the result of a capture experiment (binned into even genomic windows)
    and plots a heatmap of interactions over a specified genomic range. If a capture probe name 
    is not supplied the script will plot all probes present in the file.

    Heatmaps can also be normalised (--normalise) using either:

    - n_interactions: The number of cis interactions.
    - n_rf_n_interactions: Normalised to the number of restriction fragments making up both genomic bins
                           and by the number of cis interactions.
    - ice: `ICE normalisation <https://www.nature.com/articles/nmeth.2148>` followed by number of cis interactions
            correction.  

    \f
    Args:
        cooler_fn (os.PathLike): Path to capture cooler file containing interactions
        coordinates (Union[str, os.PathLike]): Coordinates for plotting. Either chrX:1000-2000 or bed file. 
                                               If a bed file, the capture probe name must be contained in the name.
        resolution (int): Genomic resolution to plot. Must be present in the cooler. 
        capture_names (Tuple, optional): Capture probes to plot. If None will plot all probes. Defaults to None.
        normalisation (str, optional): Normalisation for heatmap. Choose from (n_interactions|n_rf_n_interactions|ice). Defaults to None.
        cmap (str, optional): Colour map to use for heatmap. Defaults to "jet".
        vmax (float, optional): vmaxold for heatmap. Defaults to 1.
        output_prefix (os.PathLike, optional): Output prefix for heatmap. Defaults to "".
    """
    


    # Extract a bedtool object from coordinates
    bt_coords = format_coordinates(coordinates)

    # If no capture name supplied will generate plots for all in file.
    if not capture_names:
        print('Extracting probe names')
        clrs = cooler.fileops.list_coolers(cooler_fn)
        capture_names = {clr.split("/")[1] for clr in clrs}
    
    for capture in capture_names:
        for res in resolution:
            for interval in bt_coords:

                interval_name, interval_coords = convert_interval_to_coords(interval, named=True)

                if capture in interval_name:
                    
                    # Extract matrix in correct format
                    ccm = CCMatrix(cooler_fn, binsize=res, capture_name=capture, remove_capture=remove_capture)

                    if normalisation:
                        matrix = ccm.get_matrix_normalised(
                            coordinates=interval_coords, normalisation_method=normalisation
                        )
                    else:
                        matrix = ccm.get_matrix(coordinates=interval_coords)

                    print(f'Plotting {capture} at {res}')

                    fig = plot_matrix(matrix, cmap=cmap, vmax=vmax, vmin=vmin)
                    fig.savefig(f"{output_prefix}_{capture}_{res}.png")
