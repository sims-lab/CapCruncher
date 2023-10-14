# ruff: noqa: F821

import os
import sys
import pandas as pd
import numpy as np
import subprocess
import pyarrow.dataset as ds
import pathlib
from loguru import logger


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
            "--statistics",
            snakemake.output.statistics,
        ]

        with open(snakemake.log[0], "w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=f)

    else:
        logger.warning("The input dataset is empty, skipping deduplication step.")

        outdir = pathlib.Path(snakemake.output.slices)

        logger.warning(f"Creating empty output directory: {outdir}")
        outdir.mkdir(parents=True, exist_ok=True)

        logger.warning(f"Creating empty stats file: {snakemake.output.stats_read}")
        pd.DataFrame().to_csv(snakemake.output.stats)


except Exception as e:
    print(e)
    os.makedirs(snakemake.output.slices, exist_ok=True)
