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
from capcruncher.api.statistics import FastqDeduplicationStatistics
from capcruncher_tools.api import deduplicate_fastq
import pandas as pd
import pathlib



def deduplicate(
    fastq_1: List[str],
    fastq_2: List[str],
    output_prefix: Union[str, pathlib.Path] = "deduplicated_",
    statistics: str = "deduplication_statistics.json",
    sample_name: str = "sampleX",
    shuffle: bool = False,
    **kwargs,
):


    df_stats = deduplicate_fastq(
        fastq1=fastq_1,
        fastq2=fastq_2,
        output_prefix=output_prefix,
        sample_name=sample_name,
        shuffle=shuffle,
    )
    
    dedup_stats = FastqDeduplicationStatistics(
        sample=sample_name,
        total=df_stats.query("stat_type == 'read_pairs_total'")["stat"].values[0],
        duplicates=df_stats.query("stat_type == 'read_pairs_duplicated'")["stat"].values[0],
    )
    with open(statistics, "w") as f:
        f.write(dedup_stats.model_dump_json())

    

    logging.info("Printing deduplication statistics to stdout")
    # Print stats to stdout
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(tabulate.tabulate(df_vis, headers="keys", tablefmt="psql", showindex=False))
