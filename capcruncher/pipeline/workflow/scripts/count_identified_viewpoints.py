import os
import sys
import pandas as pd
import numpy as np
import subprocess

from capcruncher import api
import ibis


ibis.options.interactive = False

con = ibis.duckdb.connect(threads=snakemake.threads)
tbl = con.register(f"parquet://{snakemake.params.slices_dir}", table_name="reporters")
unique_viewpoints = (tbl[["viewpoint", "pe"]]
                        .distinct()
                        .execute(limit=None)
                        .replace("", pd.NA)
                        .dropna()
)

unique_viewpoints.to_csv(snakemake.output[0], sep="\t", index=False)


