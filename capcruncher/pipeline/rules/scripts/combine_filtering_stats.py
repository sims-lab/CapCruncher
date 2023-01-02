from capcruncher.api.statistics import collate_read_data, collate_slice_data

# Collate data
df_read_data = collate_read_data(snakemake.input.read_level_stats)
df_slice_data = collate_slice_data(snakemake.input.slice_level_stats)

# Write data
df_read_data.to_csv(snakemake.output.read_data, sep=",", index=False)
df_slice_data.to_csv(snakemake.output.slice_data, sep=",", index=False)
