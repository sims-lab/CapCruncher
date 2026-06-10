"""Genomics utilities: restriction enzymes, GTF/coordinate conversion."""

from __future__ import annotations

import pandas as pd


def get_human_readable_number_of_bp(bp: int) -> str:
    """Converts integer into human readable basepair number"""

    if bp < 1000:
        return f"{bp}bp"
    if (bp / 1e3) < 1000:
        return f"{bp / 1e3}kb"
    return f"{bp / 1e6}mb"


def convert_interval_to_coords(
    interval: dict | pd.Series, named: bool = False
) -> tuple[str, str]:
    """Converts interval object to standard genomic coordinates (chr:start-end)."""
    chrom = interval.get("chrom", interval.get("Chromosome"))
    start = interval.get("start", interval.get("Start"))
    end = interval.get("end", interval.get("End"))
    name = interval.get("name", interval.get("Name", "Unnammed"))

    if not named:
        return ("Unnammed", f"{chrom}:{start}-{end}")
    else:
        return (name, f"{chrom}:{start}-{end}")


def gtf_line_to_bed12_line(df: pd.DataFrame) -> str:
    df = df.sort_values(["seqname", "start"])
    geneid = df["geneid"].iloc[0]
    exons = df.query('feature == "exon"')
    chrom = df["seqname"].iloc[0]
    start = str(df["start"].min())
    end = str(df["end"].max())
    strand = df["strand"].iloc[0]
    thick_start = start if strand == "+" else end
    thick_end = thick_start
    color = "0,0,0"
    block_count = str(exons.shape[0])
    block_sizes = ",".join((exons["end"] - exons["start"]).values.astype(str))
    block_starts = ",".join((exons["start"] - int(start)).astype(str))

    return "\t".join(
        [
            chrom,
            start,
            end,
            geneid,
            "0",
            strand,
            thick_start,
            thick_end,
            color,
            block_count,
            block_sizes,
            block_starts,
        ]
    )


def get_restriction_site(restriction_enzyme: str) -> str:
    """Gets the restriction site for a given restriction enzyme.

    Can be either the name of the restriction enzyme or the restriction site itself.
    The restriction site will just be returned if it is a valid DNA sequence.

    Args:
        restriction_enzyme: Name of restriction enzyme or restriction site.

    Returns:
        Restriction site.

    Raises:
        ValueError: If restriction enzyme is not found.
    """
    import re

    if re.match(r"^[ACGTacgt]+$", restriction_enzyme):
        return restriction_enzyme

    import Bio.Restriction

    all_enzymes = {e.lower(): e for e in Bio.Restriction.AllEnzymes.as_string()}
    if restriction_enzyme.lower() not in all_enzymes:
        raise ValueError(f"Restriction enzyme {restriction_enzyme} not found.")
    else:
        return Bio.Restriction.AllEnzymes.get(
            all_enzymes[restriction_enzyme.lower()]
        ).site
