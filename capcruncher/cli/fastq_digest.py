import os
from typing import Literal, Tuple, Dict

import pandas as pd
from loguru import logger as logging
import polars as pl

def digest(
    fastqs: Tuple,
    restriction_site: str,
    mode: Literal["flashed", "pe"] = "pe",
    output_file: os.PathLike = "out.fastq.gz",
    minimum_slice_length: int = 18,
    statistics: os.PathLike = "",
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
        statstics: Output prefix for stats file.
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
        restriction_site=get_restriction_site(restriction_site),
        output=output_file,
        read_type=mode.title(),
        sample_name=sample_name,
        minimum_slice_length=minimum_slice_length,
    )

    logging.info("Digestion complete. Generating statistics")
    with open(statistics, "w") as f:
        f.write(stats.model_dump_json())

    return stats
