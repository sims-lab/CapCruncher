import os
from typing import Literal, Tuple, Dict

import pandas as pd
from loguru import logger as logging
import polars as pl
from capcruncher.api.statistics import FastqDigestionStatistics, HistogramLength, HistogramNumber


def digest(
    fastqs: Tuple,
    restriction_site: str,
    mode: Literal["flashed", "pe"] = "pe",
    output_file: os.PathLike = "out.fastq.gz",
    minimum_slice_length: int = 18,
    stats_prefix: os.PathLike = "",
    sample_name: str = "sampleX",
    **kwargs,
) -> Dict[str, pl.DataFrame]:
    """
    Digest FASTQ files.

    Args:
        fastqs: Tuple of FASTQ files.
        restriction_site: Restriction enzyme name or sequence to use for in silico digestion.
        mode: Digestion mode. Combined (Flashed) or non-combined (PE) read pairs.
        output_file: Output file path.
        minimum_slice_length: Minimum slice length.
        stats_prefix: Output prefix for stats file.
        sample_name: Name of sample e.g. DOX_treated_1. Required for correct statistics.

    Returns:
        A dictionary of stats: stats_read_level, stats_hist_unfilt, stats_hist_filt
    """
    from capcruncher_tools.api import digest_fastq
    from capcruncher.utils import get_restriction_site

    logging.info("Digesting FASTQ files")

    if len(fastqs) > 1 and mode == "flashed":
        raise ValueError("Flashed mode can only be used with a single FASTQ file")

    stats = digest_fastq(
        fastqs=fastqs,
        restriction_enzyme=get_restriction_site(restriction_site),
        output=output_file,
        read_type=mode.title(),
        sample_name=sample_name,
        minimum_slice_length=minimum_slice_length,
    )

    logging.info("Digestion complete")
    
    
    
    
    
    
    # histograms = dict()
    # for hist_type in ["unfilt", "filt"]:
    #     df = stats[f"stats_hist_{hist_type}"].to_pandas()
    #     df = df.rename(columns={"slice_number": "n_slices"})
    #     df = df.replace({r"^One$": 1, r"^Two$": 2}, regex=True)
    #     df = df.assign(
    #         filtered=hist_type == "filt",
    #     )
    #     histograms[hist_type] = HistogramNumber.from_dataframe(
    #         df=df,
    #         statistic_name="n_slices",
    #         read_in_pair="unpaired",
    #     )
    
    # digestion_stats = FastqDigestionStatistics(
    #     id=sample_name,
    #     read_pairs_input=
    # )

    # df_stats_read_fmt.to_csv(f"{stats_prefix}.digestion.read.summary.csv", index=False)

    return stats
