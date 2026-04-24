"""Check that there is only one restriction fragment per viewpoint."""

import pathlib
from loguru import logger

from capcruncher.api.annotate import annotate_intervals


def check_n_bins_per_viewpoint(
    *,
    bins,
    viewpoints,
    output_sentinel,
    ignore_multiple_bins_per_viewpoint,
):
    logger.info("Checking that there is only one restriction fragment per viewpoint")

    gr = annotate_intervals(
        query=viewpoints,
        annotations=bins,
        name="restriction_fragments",
        method="count",
        fraction=0.51,
    )
    multiple_fragments = (gr["restriction_fragments"] > 1).sum()
    has_multiple_frags = multiple_fragments > 0

    if has_multiple_frags and not ignore_multiple_bins_per_viewpoint:
        raise ValueError(
            f"""The following viewpoints overlap multiple restriction fragments:\n{gr}\n"""
        )

    pathlib.Path(output_sentinel).touch()


def main(snakemake):
    with logger.catch():
        check_n_bins_per_viewpoint(
            bins=snakemake.input.bins,
            viewpoints=snakemake.input.viewpoints,
            output_sentinel=snakemake.output.sentinel,
            ignore_multiple_bins_per_viewpoint=snakemake.params.ignore_multiple_bins_per_viewpoint,
        )


if "snakemake" in globals():
    main(globals()["snakemake"])
