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
    statistics_path: str = "deduplication_statistics.csv",
    sample_name: str = "sampleX",
    shuffle: bool = False,
    **kwargs,
):
    from capcruncher_tools.api import deduplicate_fastq
    import pandas as pd
    import pathlib

    print("running")

    output_prefix = pathlib.Path(output_prefix).parent
    output_prefix.mkdir(exist_ok=True, parents=True)

    stats_path = pathlib.Path(statistics_path).with_suffix(".deduplication.csv")
    stats_path.parent.mkdir(exist_ok=True, parents=True)

    if not stats_path.parent.exists():
        raise ValueError(f"Statistics path {stats_path.parent} does not exist")

    df_stats = deduplicate_fastq(
        fastq1=fastq_1,
        fastq2=fastq_2,
        output_prefix=output_prefix,
        sample_name=sample_name,
        shuffle=shuffle,
    )

    # logging.info("Saving deduplication statistics")
    # # df_stats = (
    # #     pd.Series(deduplication_results)
    # #     .to_frame("stat")
    # #     .reset_index()
    # #     .rename(columns={"index": "stat_type"})
    # #     .assign(
    # #         read_number=0,
    # #         read_type="pe",
    # #         stage="deduplication",
    # #         sample=sample_name,
    # #     )
    # )

    logging.info(f"Saving stats to {stats_path}")
    df_stats.to_csv(stats_path, index=False)

    logging.info("Printing deduplication statistics to stdout")
    # Print stats to stdout
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(tabulate.tabulate(df_vis, headers="keys", tablefmt="psql", showindex=False))
