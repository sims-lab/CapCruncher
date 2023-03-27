# ruff: noqa: F821
from capcruncher.api.statistics import collate_read_data

df = collate_read_data(snakemake.input.read_level_stats)
df.to_csv(snakemake.output[0], sep=",", index=False)
