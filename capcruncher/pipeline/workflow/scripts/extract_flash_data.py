# ruff: noqa: F821


import os
import sys
import pandas as pd
import ujson

from capcruncher.api.statistics import FlashStats

df_stats = pd.read_csv(snakemake.input[0], sep="\t")
df_stats["sample"] = df_stats["Sample"].str.split("_part").str[0]
df_stats = df_stats[["sample", "combopairs", "uncombopairs"]].groupby("sample").sum().reset_index()

stats = []
for index, row in df_stats.iterrows():
    stat = FlashStats(
        sample=row["sample"],
        n_combined=row["combopairs"],
        n_uncombined=row["uncombopairs"],)
    stats.append(stat)

with open(snakemake.output[0], "w") as f:
    stats_json = [s.model_dump_json() for s in stats]
    f.write(ujson.dumps(stats_json, indent=4))