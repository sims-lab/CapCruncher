import os

import click
import typer

plot_app = typer.Typer(
    help="Generate plots from CapCruncher outputs.",
    context_settings={"help_option_names": ["-h", "--help"]},
    invoke_without_command=True,
)


def render_plot(
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


def _render_or_raise(region: str | None, template: os.PathLike | None, output: str):
    if region is None or template is None:
        raise click.UsageError("Missing option '-r' / '--region' or '-t' / '--template'.")

    render_plot(region=region, template=template, output=output)


@plot_app.callback(invoke_without_command=True)
def plot(
    ctx: typer.Context,
    region: str | None = typer.Option(
        None,
        "-r",
        "--region",
        help="Genomic coordinates of the region to plot.",
    ),
    template: str | None = typer.Option(
        None,
        "-t",
        "--template",
        help="TOML file containing the template for the plot.",
    ),
    output: str = typer.Option(
        "capcruncher_plot.png",
        "-o",
        "--output",
        help="Output file path. The file extension determines the output format.",
    ),
):
    """Generate plots for outputs produced by CapCruncher."""

    if ctx.invoked_subcommand is not None:
        return

    _render_or_raise(region=region, template=template, output=output)


@plot_app.command("render")
def plot_render(
    region: str = typer.Option(
        ...,
        "-r",
        "--region",
        help="Genomic coordinates of the region to plot.",
    ),
    template: str = typer.Option(
        ...,
        "-t",
        "--template",
        help="TOML file containing the template for the plot.",
    ),
    output: str = typer.Option(
        "capcruncher_plot.png",
        "-o",
        "--output",
        help="Output file path. The file extension determines the output format.",
    ),
):
    """Render a PlotNado TOML template."""

    render_plot(region=region, template=template, output=output)


cli = typer.main.get_command(plot_app)
