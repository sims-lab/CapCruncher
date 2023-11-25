import click
import os
import pathlib
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import capcruncher.api.plotting as cp
from loguru import logger


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
