#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""
from typing import List, Tuple, Union
from loguru import logger as logging
import tabulate
import pathlib


def deduplicate(
    fastq_1: List[str],
    fastq_2: List[str],
    output_prefix: Union[str, pathlib.Path] = "deduplicated_",
    statistics: str = "deduplication_statistics.csv",
    sample_name: str = "sampleX",
    shuffle: bool = False,
    **kwargs,
):
    from capcruncher_tools.api import deduplicate_fastq
    import pandas as pd
    import pathlib

    df_stats = deduplicate_fastq(
        fastq1=fastq_1,
        fastq2=fastq_2,
        output_prefix=output_prefix,
        sample_name=sample_name,
        shuffle=shuffle,
    )

    logging.info(f"Saving stats to {statistics}")
    df_stats.to_csv(statistics, index=False)

    logging.info("Printing deduplication statistics to stdout")
    # Print stats to stdout
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(tabulate.tabulate(df_vis, headers="keys", tablefmt="psql", showindex=False))
