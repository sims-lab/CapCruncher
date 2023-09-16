#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""
from typing import List, Tuple
from loguru import logger as logging
import tabulate


def deduplicate(
    fastq_input: List[Tuple[str, str]],
    fastq_output: List[Tuple[str, str]],
    statistics_path: str,
    sample_name: str,
    shuffle: bool = False,
    **kwargs,
):
    from capcruncher_tools.api import deduplicate_fastq
    import pandas as pd
    import pathlib

    output_prefix = pathlib.Path(fastq_output[0][0]).parent
    output_prefix.mkdir(exist_ok=True, parents=True)

    stats_path = pathlib.Path(statistics_path).with_suffix(".deduplication.csv")
    stats_path.parent.mkdir(exist_ok=True, parents=True)

    if not stats_path.parent.exists():
        raise ValueError(f"Statistics path {stats_path.parent} does not exist")

    deduplication_results = deduplicate_fastq(fastq_input, fastq_output, shuffle)
    logging.info("Saving deduplication statistics")
    df_stats = (
        pd.Series(deduplication_results)
        .to_frame("stat")
        .reset_index()
        .rename(columns={"index": "stat_type"})
        .assign(
            read_number=0,
            read_type="pe",
            stage="deduplication",
            sample=sample_name,
        )
    )

    logging.info(f"Saving stats to {stats_path}")
    df_stats.to_csv(stats_path, index=False)

    logging.info("Printing deduplication statistics to stdout")
    # Print stats to stdout
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(tabulate.tabulate(df_vis, headers="keys", tablefmt="psql", showindex=False))
