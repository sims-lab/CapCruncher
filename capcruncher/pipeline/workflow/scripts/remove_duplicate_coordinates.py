# ruff: noqa: F821

import os
import sys
import pandas as pd
import numpy as np
import subprocess
import pyarrow.dataset as ds


# Check if the input dataset is not empty
try:
    dataset = ds.dataset(snakemake.input.slices_directory, format="parquet")
    n_rows = dataset.count_rows()

    if n_rows != 0:
        cmd = [
            "capcruncher",
            "interactions",
            "deduplicate",
            snakemake.input.slices_directory,
            "-o",
            snakemake.output.slices,
            "--read-type",
            snakemake.params.read_type,
            "--sample-name",
            snakemake.params.sample_name,
            "--stats-prefix",
            snakemake.params.stats_prefix,
        ]

        with open(snakemake.log[0], "w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=f)

    else:
        print("The input dataset is empty, skipping deduplication step.")

except Exception as e:
    print(e)
    os.makedirs(snakemake.output.slices, exist_ok=True)
