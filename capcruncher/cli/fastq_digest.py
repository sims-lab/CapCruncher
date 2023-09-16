import os
from typing import Literal, Tuple

import pandas as pd
from loguru import logger as logging


def digest(
    fastqs: Tuple,
    restriction_site: str,
    mode: Literal["flashed", "pe"] = "pe",
    output_file: os.PathLike = "out.fastq.gz",
    minimum_slice_length: int = 18,
    stats_prefix: os.PathLike = "",
    sample_name: str = "sampleX",
    **kwargs,
):
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
        A dictionary of stats.
    """
    from capcruncher_tools.api import digest_fastq

    logging.info("Digesting FASTQ files")

    stats = digest_fastq(
        fastqs=fastqs,
        restriction_enzyme=restriction_site,
        output=output_file,
        read_type=mode.title(),
        sample_name=sample_name,
        minimum_slice_length=minimum_slice_length,
    )

    logging.info("Digestion complete")
    df_stats_read = stats["stats_read_level"].to_pandas()
    df_stats_read_fmt = pd.DataFrame(
        {
            "stat_type": ["unfiltered", "filtered"],
            "stat": [
                df_stats_read["number_of_read_pairs_unfiltered"].sum(),
                df_stats_read["number_of_read_pairs_filtered"].sum(),
            ],
        }
    ).assign(stage="digestion", sample=sample_name, read_type=mode)

    for hist_type in ["unfilt", "filt"]:
        df = stats[f"stats_hist_{hist_type}"].to_pandas()
        df = df.rename(columns={"slice_number": "n_slices"})
        df = df.replace({r"^One$": 1, r"^Two$": 2}, regex=True)
        df.to_csv(f"{stats_prefix}.digestion.{hist_type}.histogram.csv", index=False)

    df_stats_read_fmt.to_csv(f"{stats_prefix}.digestion.read.summary.csv ", index=False)

    return stats
