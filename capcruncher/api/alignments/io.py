from __future__ import annotations

from collections import namedtuple
import os
from pathlib import Path

from loguru import logger
import pandas as pd
import pysam
import xxhash

CCAlignment = namedtuple(
    "CCAlignment",
    field_names=[
        "slice_id",
        "slice_name",
        "parent_id",
        "parent_read",
        "pe",
        "slice",
        "uid",
        "mapped",
        "multimapped",
        "chrom",
        "start",
        "end",
        "coordinates",
    ],
)


def parse_alignment(aln: pysam.AlignedSegment) -> CCAlignment:
    """Parses reads from a bam file into a list.

    Extracts:
     -read name
     -parent reads
     -flashed status
     -slice number
     -mapped status
     -multimapping status
     -chromosome number (e.g. chr10)
     -start (e.g. 1000)
     -end (e.g. 2000)
     -coords e.g. (chr10:1000-2000)


    Args:
     aln: pysam.AlignmentFile.
    Returns:
     list: Containing the attributes extracted.

    """

    slice_name = aln.query_name
    parent_read, pe, slice_number, uid = slice_name.split("|")
    parent_id = xxhash.xxh3_64_intdigest(parent_read, seed=42)
    slice_id = xxhash.xxh3_64_intdigest(slice_name, seed=42)
    ref_name = aln.reference_name
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    # Check if read mapped
    if aln.is_unmapped:
        mapped = 0
        multimapped = 0
        ref_name = ""
        ref_start = 0
        ref_end = 0
        coords = ""
    else:
        mapped = 1
        coords = f"{ref_name}:{ref_start}-{ref_end}"
        # Check if multimapped
        if aln.is_secondary:
            multimapped = 1
        else:
            multimapped = 0

    return CCAlignment(
        slice_id=slice_id,
        slice_name=slice_name,
        parent_id=parent_id,
        parent_read=parent_read,
        pe=pe.lower(),
        slice=int(slice_number),
        uid=int(uid),
        mapped=mapped,
        multimapped=multimapped,
        chrom=ref_name,
        start=int(ref_start),
        end=int(ref_end),
        coordinates=coords,
    )


def parse_bam(bam: Path | str) -> pd.DataFrame:
    """Uses parse_alignment function convert bam file to a dataframe.

    Extracts:
     -'slice_name'
     -'parent_read'
     -'pe'
     -'slice'
     -'mapped'
     -'multimapped'
     -'chrom'
     -'start'
     -'end'
     -'coordinates'

    Args:
        bam: Path to bam file.

    Returns:
     pd.Dataframe: DataFrame with the columns listed above.

    """

    import numpy as np

    bam = os.fspath(bam)

    # Load reads into dataframe
    logger.info("Parsing BAM file")
    df_bam = pd.DataFrame(
        [
            parse_alignment(aln)
            for aln in pysam.AlignmentFile(bam, "rb").fetch(until_eof=True)
        ],
    )
    df_bam["bam"] = os.path.basename(bam)

    # Perform dtype conversions
    logger.info("Converting dtypes")
    df_bam["chrom"] = df_bam["chrom"].astype("category")
    pe_category = pd.CategoricalDtype(["flashed", "pe"])
    df_bam["pe"] = df_bam["pe"].astype(
        pe_category
    )  # Only the one type present so need to include both
    df_bam["coordinates"] = df_bam["coordinates"].astype("category")
    df_bam["parent_read"] = df_bam["parent_read"].astype("category")
    df_bam["slice"] = df_bam["slice"].astype(np.int8)
    df_bam["uid"] = df_bam["uid"].astype(np.int8)
    df_bam["multimapped"] = df_bam["multimapped"].astype(bool)
    df_bam["mapped"] = df_bam["mapped"].astype(bool)
    df_bam["bam"] = df_bam["bam"].astype("category")

    logger.info("Finished parsing BAM file")
    return df_bam


def bam_to_parquet(
    bam: Path | str, output: Path | str
) -> Path | str:
    """Converts bam file to parquet file.

    Args:
     bam: Path to bam file.
     output: Path to output parquet file.

    """
    df_bam = parse_bam(bam)
    df_bam.to_parquet(output)

    return output
