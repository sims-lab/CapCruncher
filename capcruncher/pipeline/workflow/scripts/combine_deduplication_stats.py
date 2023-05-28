# ruff: noqa: F821
from capcruncher.api.statistics import collate_read_data

# Read in the data
df = collate_read_data(snakemake.input)

# # Ignore the reads_removed stat_type
# df = df[df["stat_type"] != "reads_removed"]

# Write out the data
df.to_csv(snakemake.output[0], sep=",", index=False)
