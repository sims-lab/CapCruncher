import os


def plot(
    region: str,
    template: os.PathLike,
    output: str,
) -> None:
    """Plot a region using a template.

    Args:
        region (str): Genomic region to plot.
        template (os.PathLike): Path to template file.
        output (str): Path to output file.

    """
    from plotnado import GenomicFigure

    fig = GenomicFigure.from_toml(str(template))
    fig.save(output, region=region)
