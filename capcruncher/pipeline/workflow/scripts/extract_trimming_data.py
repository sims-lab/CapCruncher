# ruff: noqa: F821


import os
import sys
import pandas as pd
import ujson

from capcruncher.api.statistics import FastqTrimmingStatistics


df_stats = pd.read_csv(snakemake.input[0], sep="\t")
df_stats["read_number"] = df_stats["Sample"].str.split("_").str[-1].astype(int)
df_stats["sample"] = df_stats["Sample"].str.extract(r"(.+)_part\d+_\d+").iloc[:, 0]
df_stats_agg = df_stats.groupby(["sample", "read_number"]).sum().reset_index()

stats = []
for index, row in df_stats_agg.iterrows():
    stat = FastqTrimmingStatistics.from_multiqc_entry(row)
    stats.append(stat)
    

with open(snakemake.output[0], "w") as f:
    stats_json = [s.model_dump_json() for s in stats]
    f.write(ujson.dumps(stats_json, indent=4))