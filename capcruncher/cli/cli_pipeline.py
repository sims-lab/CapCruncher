import os
from capcruncher.cli import cli
import click
from importlib import import_module, metadata
import subprocess
import sys
import logging

@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option("-h", "--help", "show_help", is_flag=True)
@click.option("--version", "show_version", is_flag=True)
@click.version_option(metadata.version(distribution_name="capcruncher"))
@click.argument(
    "mode",
    type=click.Choice(["make", "run", "plot", "show", "clone", "touch", "report"]),
)
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
@click.option("-s", "--pipeline-statistics-path", help="Path to capcruncher_statistics directory", default='capcruncher_statistics')
@click.option("-o", "--pipeline-report-path", help="Output path for CapCruncher report", default='capcruncher_report.html')
def pipeline(mode, pipeline_options, show_help=False, show_version=False, **report_options):

    """Runs the data processing pipeline"""

    fn = os.path.abspath(__file__)
    dir_cli = os.path.dirname(fn)
    dir_package = os.path.dirname(dir_cli)

    if mode in ["make", "run", "plot", "show", "clone", "touch"]:

        cmd = [
            "python",
            f"{dir_package}/pipeline/pipeline.py",
        ]

        if show_help:
            cmd.append("--help")
            subprocess.run(cmd)
            sys.exit()

        cmd.append(mode.replace("run", "make"))

        if pipeline_options:
            cmd.extend(pipeline_options)

        # Implicitly deal with the missing --local option
        if (
            not os.path.exists(os.environ.get("DRMAA_LIBRARY_PATH", ""))
            and not "--local" in pipeline_options
        ):
            logging.warning(
                "DRMAA_LIBRARY_PATH is incorrect. Implicitly using --local with 4 cores"
            )
            cmd.append("--local")
            cmd.append("-p 4")

        completed = subprocess.run(cmd)

        if not completed.returncode == 0:
            raise RuntimeError(
                "CapCruncher pipeline failed. Check pipeline.log for details"
            )

    else:
        from capcruncher.pipeline.make_report import generate_report
        generate_report(**report_options)
        