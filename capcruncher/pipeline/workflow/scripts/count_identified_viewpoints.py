import os
import sys
import pandas as pd
import numpy as np
import subprocess

from capcruncher import api
import ibis


ibis.options.interactive = False

con = ibis.duckdb.connect(threads=snakemake.threads)
parquet_path = snakemake.params.slices_dir
if os.path.isdir(parquet_path):
    parquet_path = os.path.join(parquet_path, "*.parquet")

tbl = con.read_parquet(parquet_path)
unique_viewpoints = (tbl[["viewpoint", "pe"]]
                        .distinct()
                        .execute()
                        .replace("", pd.NA)
                        .dropna()
)

unique_viewpoints.to_csv(snakemake.output[0], sep="\t", index=False)

