from pathlib import Path
from typing import Any

from loguru import logger


def digest_genome(
    input_fasta: Path | str,
    recognition_site: str,
    output_file: Path | str = "genome_digest.bed",
    sort: bool = False,
    remove_cutsite: bool = True,
    **kwargs: Any,
) -> None:
    """Digest a genome FASTA and optionally sort the resulting BED file."""

    import polars as pl
    from capcruncher_tools.api import digest_genome as digest_genome_records

    from capcruncher.utils import get_restriction_site

    logger.info("Digesting genome")
    digest_genome_records(
        fasta=str(input_fasta),
        output=str(output_file),
        restriction_enzyme=get_restriction_site(recognition_site),
        remove_recognition_site=remove_cutsite,
        minimum_slice_length=18,
        n_threads=1,
    )

    logger.info("Digestion complete")

    if sort:
        logger.info("Sorting output")
        df = pl.read_csv(
            output_file,
            separator="\t",
            new_columns=["chrom", "start", "end", "name"],
            schema_overrides={
                "chrom": pl.String,
                "start": pl.Int64,
                "end": pl.Int64,
                "name": pl.String,
            },
        )

        df = (
            df.filter(pl.col("end") > pl.col("start"))
            .sort(["chrom", "start"])
            .drop(["name"])
            .with_row_index("name")[["chrom", "start", "end", "name"]]
        )

        df.write_csv(output_file, separator="\t", include_header=False)
