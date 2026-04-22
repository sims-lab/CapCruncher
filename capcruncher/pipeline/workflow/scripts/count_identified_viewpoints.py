import os
import sys
import duckdb
import pandas as pd
import numpy as np
import subprocess

from capcruncher import api


con = duckdb.connect()
con.execute(f"PRAGMA threads={int(snakemake.threads)}")
parquet_path = str(snakemake.params.slices_dir).replace("'", "''")
unique_viewpoints = con.sql(
        f"""
        SELECT DISTINCT viewpoint, pe
        FROM read_parquet('{parquet_path}')
        WHERE nullif(viewpoint, '') IS NOT NULL
            AND nullif(pe, '') IS NOT NULL
        """
).df()

unique_viewpoints.to_csv(snakemake.output[0], sep="\t", index=False)


