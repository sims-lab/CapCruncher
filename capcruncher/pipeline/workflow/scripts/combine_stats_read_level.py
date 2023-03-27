# ruff: noqa: F821
import pandas as pd

from capcruncher.utils import read_dataframes


dframes = read_dataframes(snakemake.input)
df = pd.concat(dframes)
df.sort_values(["sample", "read_type", "stat"], ascending=[True, True, False]).to_csv(
    snakemake.output[0], sep=",", index=False
)
