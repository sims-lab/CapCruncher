import os
from capcruncher.cli import cli
import click
from importlib import metadata
import subprocess
import sys


@cli.command(context_settings=dict(ignore_unknown_options=True), name="pipeline")
@click.option("-h", "--help", "show_help", is_flag=True)
@click.option("--version", "show_version", is_flag=True)
@click.version_option(metadata.version(distribution_name="capcruncher"))
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def pipeline(pipeline_options, show_help=False, show_version=False):

    """Runs the data processing pipeline"""

    fn = os.path.abspath(__file__)
    dir_cli = os.path.dirname(fn)
    dir_package = os.path.dirname(dir_cli)

    cmd = [
        "snakemake",
        "-s",
        f"{dir_package}/pipeline/workflow/Snakefile",
    ]

    if show_help:
        # Run snakemake with --help
        # Capture the output and replace usage: snakemake with usage: capcruncher pipeline
        # Print the output
        cmd.append("--help")
        _completed = subprocess.run(cmd, capture_output=True, shell=False)
        output = _completed.stdout.decode("utf-8")
        output = output.replace("usage: snakemake", "usage: capcruncher pipeline")
        click.echo(f"\n{output}")
        sys.exit(0)

    if pipeline_options:
        excluded_options = ["--version", "make", "run", "show"]

        cmd.extend(
            [option for option in pipeline_options if option not in excluded_options]
        )

    # Implicitly deal with a missing --cores option
    if "--cores" not in pipeline_options and "-c" not in pipeline_options:
        cmd.append("--cores 1")

    _completed = subprocess.run(cmd)


@cli.command(name="pipeline-config")
@click.option("-h", "--help", "show_help", is_flag=True)
@click.option("--version", "show_version", is_flag=True)
@click.version_option(metadata.version(distribution_name="capcruncher"))
@click.option(
    "-i", "--input", "input_files", type=click.Path(exists=True), multiple=True
)
@click.option("--generate-design", is_flag=True)
def pipeline_config(*args, **kwargs):
    """Configures the data processing pipeline"""

    from cookiecutter.main import cookiecutter
    import pathlib

    fn = pathlib.Path(__file__).resolve()
    dir_cli = fn.parent
    dir_package = dir_cli.parent

    cookiecutter(str(dir_package / "pipeline" / "config"))
