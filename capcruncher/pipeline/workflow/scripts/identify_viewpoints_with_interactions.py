# ruff: noqa: F821
import os
import sys
import cooler
import pandas as pd
from collections import defaultdict
import ujson as json
import pathlib


coolers = snakemake.input[0]
samples = snakemake.params.samples

for (sample, clr) in zip(samples, coolers):
    viewpoints_per_sample = []
    # Check all groups in the cooler file to see if they have any counts
    viewpoints = cooler.api.list_coolers(clr)

    for viewpoint in viewpoints:
        clr = cooler.Cooler(f"{clr}::{viewpoint}")
        count = clr.pixels()[:].sum()

        if count > 0:
            viewpoints_per_sample.append(viewpoint)

    # Save the viewpoints with counts to a file
    with open(pathlib.Path(snakemake.params.outdir) / f"{sample}.json", "w") as f:
        json.dump(viewpoints_per_sample, f)
