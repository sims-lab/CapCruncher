# ruff: noqa: F821
import pandas as pd
from capcruncher.api.statistics import collate_read_data, collate_histogram_data

# Collate data
df_read_data = collate_read_data(snakemake.input.read_level_stats)
df_histogram_unfiltered = collate_histogram_data(snakemake.input.histogram_unfiltered)
df_histogram_filtered = collate_histogram_data(snakemake.input.histogram_filtered)

# Merge filtered and unfiltered histograms
df_hist = pd.concat(
    [
        df_histogram_unfiltered.assign(filtered=0),
        df_histogram_filtered.assign(filtered=1),
    ]
).sort_values(["sample", "read_type", "n_slices"])

# Output histogram, slice and read statics
df_hist.to_csv(snakemake.output.histogram, index=False)

# Output read data
df_read_data.to_csv(snakemake.output.read_data, sep=",", index=False)
