import pandas as pd

from capcruncher.utils import read_dataframes
from capcruncher.tools.statistics import collate_read_data

df = collate_read_data(snakemake.input.read_level_stats)
df.to_csv(snakemake.output[0], sep=',', index=False)
