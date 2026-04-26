import json

import pandas as pd

from capcruncher.api.statistics import FlashStats


def extract_flash_stats(flash_summary_path):
    df_stats = pd.read_csv(flash_summary_path, sep="\t")
    df_stats["sample"] = df_stats["Sample"].str.split("_part").str[0]
    df_stats = (
        df_stats[["sample", "combopairs", "uncombopairs"]]
        .groupby("sample")
        .sum()
        .reset_index()
    )

    return [
        FlashStats(
            sample=row["sample"],
            n_combined=row["combopairs"],
            n_uncombined=row["uncombopairs"],
        )
        for _, row in df_stats.iterrows()
    ]


def write_flash_stats(flash_summary_path, output_path):
    stats_json = [
        stat.model_dump_json() for stat in extract_flash_stats(flash_summary_path)
    ]
    with open(output_path, "w") as f:
        f.write(json.dumps(stats_json, indent=4))


if "snakemake" in globals():
    write_flash_stats(
        globals()["snakemake"].input[0],
        globals()["snakemake"].output[0],
    )
