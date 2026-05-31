import pathlib

import pandas as pd
from loguru import logger


def can_group_tracks_by_condition(design: pd.DataFrame) -> bool:
    return design.groupby("condition").size().max() > 1


def _track_table(paths, viewpoint, design):
    return (
        pd.DataFrame([pathlib.Path(path) for path in paths], columns=["fn"])
        .assign(
            samplename_and_vp=lambda df: df.fn.apply(lambda path: path.stem),
            samplename=lambda df: df.samplename_and_vp.str.replace(
                f"_{viewpoint}", "", regex=False
            ),
        )
        .merge(design, left_on="samplename", right_on="sample", how="left")
    )


def add_bigwig_tracks(fig, bigwigs, *, viewpoint, design, can_group_tracks):
    if not bigwigs:
        return

    logger.info("Adding bigwig tracks")
    if can_group_tracks:
        df_bw = _track_table(bigwigs, viewpoint, design)
        for condition, df in df_bw.groupby("condition"):
            fig.bigwig_collection(
                [str(fn) for fn in df.fn.tolist()],
                title=condition,
            )
            logger.info(f"Added {condition} bigwig track")
            fig.spacer()
        return

    for bigwig in bigwigs:
        bigwig_path = pathlib.Path(bigwig)
        fig.bigwig(
            bigwig,
            title=bigwig_path.stem,
            min_value=0,
            max_value=None,
        )
        logger.info(f"Added {bigwig_path.stem} bigwig track")
        fig.spacer()


def add_subtraction_tracks(fig, subtractions):
    if not subtractions:
        return

    logger.info("Adding subtraction tracks")
    for subtraction in subtractions:
        subtraction_path = pathlib.Path(subtraction)
        logger.info(f"Adding {subtraction_path.stem} subtraction track")
        fig.bigwig(subtraction, title=subtraction_path.stem)
        fig.spacer()


def add_heatmap_tracks(
    fig,
    heatmaps,
    *,
    viewpoint,
    design,
    can_group_tracks,
    binsize,
    normalization_method,
):
    if not heatmaps:
        return

    logger.info("Adding heatmaps")
    if can_group_tracks:
        df_hm = _track_table(heatmaps, viewpoint, design)

        for condition, df in df_hm.groupby("condition"):
            logger.info(f"Adding {condition} heatmap track")
            for heatmap in df.fn.tolist():
                fig.add_track(
                    "capcruncher",
                    file=str(heatmap),
                    title=f"{condition}: {pathlib.Path(heatmap).stem}",
                    resolution=binsize,
                    viewpoint=viewpoint,
                    normalisation=normalization_method,
                    balance=False,
                )
            fig.spacer()
        return

    for heatmap in heatmaps:
        heatmap_path = pathlib.Path(heatmap)
        logger.info(f"Adding {heatmap_path.stem} heatmap track")
        fig.add_track(
            "capcruncher",
            file=heatmap,
            title=heatmap_path.stem,
            resolution=binsize,
            viewpoint=viewpoint,
            normalisation=normalization_method,
            balance=False,
        )
        fig.spacer()


def build_figure(
    *,
    bigwigs,
    subtractions,
    heatmaps,
    genes,
    design,
    viewpoint,
    binsize,
    normalization_method,
):
    from plotnado import GenomicFigure

    logger.info("Checking if we can group tracks by condition")
    can_group_tracks = can_group_tracks_by_condition(design)

    logger.info("Setting up tracks")
    fig = GenomicFigure(theme=None)
    fig.scalebar()

    add_bigwig_tracks(
        fig,
        bigwigs,
        viewpoint=viewpoint,
        design=design,
        can_group_tracks=can_group_tracks,
    )
    add_subtraction_tracks(fig, subtractions)
    add_heatmap_tracks(
        fig,
        heatmaps,
        viewpoint=viewpoint,
        design=design,
        can_group_tracks=can_group_tracks,
        binsize=binsize,
        normalization_method=normalization_method,
    )

    if genes and pathlib.Path(genes).is_file():
        logger.info("Adding genes track")
        fig.genes(data=genes)
        fig.spacer()

    fig.axis()
    return fig


def save_figure(fig, *, output_fig, output_template, coordinates):
    logger.info("Making figure")

    logger.info(f"Saving figure to: {output_fig}")
    fig.save(output_fig, region=coordinates)

    logger.info(f"Exporting template to {output_template}")
    fig.to_toml(output_template)


def main(snakemake):
    logger.add(snakemake.log[0], format="{time} {level} {message}", level="INFO")

    with logger.catch():
        fig = build_figure(
            bigwigs=snakemake.input.bigwigs,
            subtractions=snakemake.input.subtractions,
            heatmaps=snakemake.input.heatmaps,
            genes=snakemake.params.genes,
            design=snakemake.params.design,
            viewpoint=snakemake.params.viewpoint,
            binsize=snakemake.params.binsize,
            normalization_method=snakemake.params.normalization_method,
        )
        save_figure(
            fig,
            output_fig=snakemake.output.fig,
            output_template=snakemake.output.template,
            coordinates=snakemake.params.coordinates,
        )


if "snakemake" in globals():
    main(globals()["snakemake"])
