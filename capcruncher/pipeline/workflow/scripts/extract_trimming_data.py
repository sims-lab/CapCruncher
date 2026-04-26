import json

import pandas as pd

from capcruncher.api.statistics import FastqTrimmingStatistics


def extract_trimming_stats(trimming_summary_path):
    df_stats = pd.read_csv(trimming_summary_path, sep="\t")
    df_stats["read_number"] = df_stats["Sample"].str.split("_").str[-1].astype(int)
    df_stats["sample"] = df_stats["Sample"].str.extract(r"(.+)_part\d+_\d+").iloc[:, 0]
    df_stats_agg = df_stats.groupby(["sample", "read_number"]).sum().reset_index()

    return [
        FastqTrimmingStatistics.from_multiqc_entry(row)
        for _, row in df_stats_agg.iterrows()
    ]


def write_trimming_stats(trimming_summary_path, output_path):
    stats_json = [
        stat.model_dump_json() for stat in extract_trimming_stats(trimming_summary_path)
    ]
    with open(output_path, "w") as f:
        f.write(json.dumps(stats_json, indent=4))


if "snakemake" in globals():
    write_trimming_stats(
        globals()["snakemake"].input[0],
        globals()["snakemake"].output[0],
    )
