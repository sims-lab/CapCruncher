import pandas as pd

from capcruncher.utils import read_dataframes
from capcruncher.tools.statistics import collate_cis_trans_data

# Collate data
df = collate_cis_trans_data(snakemake.input.cis_and_trans_stats)

# Write data
df.to_csv(snakemake.output[0], sep=",", index=False)


