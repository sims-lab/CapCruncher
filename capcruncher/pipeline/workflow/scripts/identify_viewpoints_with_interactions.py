# ruff: noqa: F821
import os
import sys
import cooler
import pandas as pd
from collections import defaultdict
import ujson as json


cooler_uri = snakemake.input[0]
viewpoints_total_counts = defaultdict(int)

# Check all groups in the cooler file to see if they have any counts
viewpoints = cooler.api.list_coolers(cooler_uri)

for viewpoint in viewpoints:
    clr = cooler.Cooler(f"{cooler_uri}::{viewpoint}")
    viewpoints_total_counts[viewpoint] = clr.pixels()[:].sum()

    if viewpoints_total_counts[viewpoint] > 0:
        # Output a dummy file to indicate that the viewpoint is present
        with open(
            f"{snakemake.params.outdir}/{snakemake.params.sample}_{viewpoint}", "w"
        ) as f:
            pass
