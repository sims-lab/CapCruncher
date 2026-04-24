# ruff: noqa: F821

import pathlib

import pandas as pd
from loguru import logger
from plotnado import GenomicFigure

logger.add(open(snakemake.log[0], "w"))


with logger.catch():
    logger.info("Checking if we can group tracks by condition")
    can_group_tracks = (
        True if snakemake.params.design.groupby("condition").size().max() > 1 else False
    )

    logger.info("Setting up tracks")
    fig = GenomicFigure()

    # Add scale bar
    fig.scalebar()

    # Bigwig tracks
    if snakemake.input.bigwigs:
        logger.info("Adding bigwig tracks")
        if can_group_tracks:
            df_bw = pd.DataFrame(
                [pathlib.Path(p) for p in snakemake.input.bigwigs], columns=["fn"]
            )
            df_bw = df_bw.assign(
                samplename_and_vp=lambda df: df.fn.apply(lambda x: x.stem),
                samplename=lambda df: df.samplename_and_vp.str.replace(
                    f"_{snakemake.params.viewpoint}", ""
                ),
            ).merge(
                snakemake.params.design,
                left_on="samplename",
                right_on="sample",
                how="left",
            )

            for condition, df in df_bw.groupby("condition"):
                fig.bigwig_collection(
                    [str(fn) for fn in df.fn.tolist()],
                    title=condition,
                )
                logger.info(f"Added {condition} bigwig track")
                fig.spacer()

        else:
            for bw in snakemake.input.bigwigs:
                bw_path = pathlib.Path(bw)
                fig.bigwig(
                    bw,
                    title=bw_path.stem,
                    min_value=0,
                    max_value=None,
                )
                logger.info(f"Added {bw_path.stem} bigwig track")
                fig.spacer()

    # Add subtractions if available
    if snakemake.input.subtractions:
        logger.info("Adding subtraction tracks")
        for sub in snakemake.input.subtractions:
            sub_path = pathlib.Path(sub)
            logger.info(f"Adding {sub_path.stem} subtraction track")
            fig.bigwig(sub, title=sub_path.stem)
            fig.spacer()

    # Add heatmaps if available
    if snakemake.input.heatmaps:
        logger.info("Adding heatmaps")
        if can_group_tracks:
            df_hm = pd.DataFrame(
                [pathlib.Path(p) for p in snakemake.input.heatmaps], columns=["fn"]
            )
            df_hm = df_hm.assign(
                samplename_and_vp=lambda df: df.fn.apply(lambda x: x.stem),
                samplename=lambda df: df.samplename_and_vp.str.replace(
                    f"_{snakemake.params.viewpoint}", ""
                ),
            ).merge(
                snakemake.params.design,
                left_on="samplename",
                right_on="sample",
                how="left",
            )

            for condition, df in df_hm.groupby("condition"):
                logger.info(f"Adding {condition} heatmap track")
                for heatmap in df.fn.tolist():
                    fig.add_track(
                        "capcruncher",
                        file=str(heatmap),
                        title=f"{condition}: {pathlib.Path(heatmap).stem}",
                        resolution=snakemake.params.binsize,
                        viewpoint=snakemake.params.viewpoint,
                        normalisation=snakemake.params.normalization_method,
                        balance=False,
                    )
                fig.spacer()
        else:
            for hm in snakemake.input.heatmaps:
                hm_path = pathlib.Path(hm)
                logger.info(f"Adding {hm_path.stem} heatmap track")
                fig.add_track(
                    "capcruncher",
                    file=hm,
                    title=hm_path.stem,
                    resolution=snakemake.params.binsize,
                    viewpoint=snakemake.params.viewpoint,
                    normalisation=snakemake.params.normalization_method,
                    balance=False,
                )
                fig.spacer()

    # Add genes if available
    if snakemake.params.genes:
        logger.info("Adding genes track")
        genes = snakemake.params.genes
        fig.genes(data=genes)
        fig.spacer()

    # Add X-axis
    fig.axis()

    # Make figure and save
    logger.info("Making figure")

    logger.info(f"Saving figure to: {snakemake.output.fig}")
    fig.save(snakemake.output.fig, region=snakemake.params.coordinates)

    # Export template used to make figure
    logger.info(f"Exporting template to {snakemake.output.template}")
    fig.to_toml(snakemake.output.template)
