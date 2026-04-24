"""
Aim: Check that there is only one restriction fragment per viewpoint.
"""

import pathlib

from loguru import logger

from capcruncher.api.annotate import annotate_intervals


with logger.catch():
    logger.info("Checking that there is only one restriction fragment per viewpoint")

    bins = snakemake.input.bins
    viewpoints = snakemake.input.viewpoints

    gr = annotate_intervals(
        query=viewpoints,
        annotations=bins,
        name="restriction_fragments",
        method="count",
        fraction=0.51,
    )
    multiple_fragments = (gr["restriction_fragments"] > 1).sum()
    has_multiple_frags = multiple_fragments > 0

    if (
        has_multiple_frags
        and not snakemake.params.ignore_multiple_bins_per_viewpoint
    ):
        raise ValueError(
            f"""The following viewpoints overlap multiple restriction fragments:\n{gr}\n"""
        )
    else:
        pathlib.Path(snakemake.output.sentinel).touch()
