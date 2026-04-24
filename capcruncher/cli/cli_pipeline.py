import os
from capcruncher.cli import cli, get_capcruncher_version
import click
from importlib import resources
import subprocess
import sys
import pathlib
import shutil


BUILTIN_PIPELINE_PRESETS = (
    "local",
    "local-conda",
    "local-apptainer",
    "slurm",
    "slurm-apptainer",
)


def has_snakemake_option(options, long_name, short_name=None):
    option_names = [long_name]
    if short_name:
        option_names.append(short_name)

    return any(
        option == option_name or option.startswith(f"{option_name}=")
        for option in options
        for option_name in option_names
    )


def get_capcruncher_config_dir() -> pathlib.Path:
    xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
    if xdg_config_home:
        return pathlib.Path(xdg_config_home).expanduser() / "capcruncher"

    return pathlib.Path.home() / ".config" / "capcruncher"


def get_pipeline_preset_dir() -> pathlib.Path:
    return get_capcruncher_config_dir() / "profiles"


def resolve_pipeline_preset(preset: str) -> pathlib.Path:
    preset_path = pathlib.Path(preset).expanduser()
    if preset_path.exists():
        return preset_path.resolve()

    bundled_path = get_pipeline_preset_dir() / preset
    if bundled_path.exists():
        return bundled_path.resolve()

    raise click.ClickException(
        f"Unknown pipeline preset '{preset}'. Run 'capcruncher pipeline-init' to install presets or pass a profile path."
    )


def install_pipeline_preset(
    preset_name: str, output_dir: pathlib.Path, force: bool
) -> pathlib.Path:
    source_dir = resources.files("capcruncher").joinpath(
        "pipeline", "profiles", preset_name
    )
    destination_dir = output_dir / preset_name

    if destination_dir.exists() and not force:
        raise click.ClickException(
            f"Preset '{preset_name}' already exists at {destination_dir}. Use --force to overwrite it."
        )

    if destination_dir.exists() and force:
        shutil.rmtree(destination_dir)

    destination_dir.mkdir(parents=True, exist_ok=True)
    for child in source_dir.iterdir():
        if child.is_file():
            (destination_dir / child.name).write_text(
                child.read_text(), encoding="utf-8"
            )

    return destination_dir


@cli.command(context_settings=dict(ignore_unknown_options=True), name="pipeline")
@click.option("-h", "--help", "show_help", is_flag=True)
@click.option(
    "--logo/--no-logo",
    default=True,
    help="Show the capcruncher logo",
    show_default=True,
)
@click.option(
    "--preset",
    type=str,
    help="CapCruncher-managed execution preset name or a profile directory path.",
)
@click.option(
    "--scale-resources",
    type=float,
    default=None,
    help="Scale workflow memory and runtime requests for retries and constrained clusters.",
)
@click.version_option(get_capcruncher_version())
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def pipeline(
    pipeline_options,
    show_help=False,
    logo=True,
    preset=None,
    scale_resources=None,
):
    """Runs the data processing pipeline"""

    fn = pathlib.Path(__file__).resolve()
    dir_cli = fn.parent
    dir_package = dir_cli.parent

    cmd = [
        "snakemake",
        "-s",
        str(dir_package / "pipeline/workflow/Snakefile"),
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

    if preset:
        if has_snakemake_option(pipeline_options, "--profile"):
            raise click.ClickException("Use either --preset or --profile, not both.")
        cmd.extend(["--profile", str(resolve_pipeline_preset(preset))])

    # Implicitly deal with a missing --cores option
    if not has_snakemake_option(pipeline_options, "--cores", "-c"):
        cmd.extend(["--cores", "1"])

    # Add the --show-failed-logs option if it is not already present
    if not has_snakemake_option(pipeline_options, "--show-failed-logs"):
        cmd.append("--show-failed-logs")

    if logo:
        with open(dir_package / "data" / "logo.txt", "r") as f:
            click.echo(f.read())

    env = os.environ.copy()
    if scale_resources is not None:
        env["SCALE_RESOURCES"] = str(scale_resources)

    # Run the pipeline
    _completed = subprocess.run(cmd, env=env)

    # If the pipeline fails, exit with the return code
    if _completed.returncode != 0:
        sys.exit(_completed.returncode)
    else:
        # Touch all files to correct timestamps
        subprocess.run(
            [
                "snakemake",
                "-s",
                str(dir_package / "pipeline/workflow/Snakefile"),
                "--touch",
                "--cores",
                "1",
            ],
            env=env,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


@cli.command(name="pipeline-init")
@click.option(
    "--output-dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=pathlib.Path),
    default=None,
    help="Directory where CapCruncher-managed pipeline presets should be installed.",
)
@click.option(
    "--preset",
    "preset_names",
    type=click.Choice(BUILTIN_PIPELINE_PRESETS),
    multiple=True,
    help="Install only the selected preset. Repeat to install multiple presets.",
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite existing preset directories if they already exist.",
)
def pipeline_init(output_dir=None, preset_names=(), force=False):
    """Installs CapCruncher-managed Snakemake presets."""

    destination_root = output_dir or get_pipeline_preset_dir()
    destination_root.mkdir(parents=True, exist_ok=True)

    presets_to_install = preset_names or BUILTIN_PIPELINE_PRESETS
    installed = []
    for preset_name in presets_to_install:
        installed.append(install_pipeline_preset(preset_name, destination_root, force))

    click.echo(f"Installed {len(installed)} pipeline preset(s) to {destination_root}")
    for installed_preset in installed:
        click.echo(f"- {installed_preset.name}: {installed_preset}")


@cli.command(name="pipeline-config")
@click.option("-h", "--help", "show_help", is_flag=True)
@click.version_option(get_capcruncher_version())
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
