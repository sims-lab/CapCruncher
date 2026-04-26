import os
from capcruncher.cli import CONTEXT_SETTINGS, cli, get_capcruncher_version
import click
from importlib import resources
import subprocess
import sys
import pathlib
import shutil
import typer


PIPELINE_PRESET_SOURCES = {
    "capcruncher-local": "local",
    "capcruncher-local-conda": "local-conda",
    "capcruncher-local-apptainer": "local-apptainer",
    "capcruncher-slurm": "slurm",
    "capcruncher-slurm-apptainer": "slurm-apptainer",
}
LEGACY_PIPELINE_PRESET_ALIASES = {
    source_name: preset_name
    for preset_name, source_name in PIPELINE_PRESET_SOURCES.items()
}
BUILTIN_PIPELINE_PRESETS = tuple(PIPELINE_PRESET_SOURCES)
PIPELINE_PRESET_CHOICES = (
    *BUILTIN_PIPELINE_PRESETS,
    *LEGACY_PIPELINE_PRESET_ALIASES,
)
PIPELINE_SUBCOMMANDS = {"run", "init", "config"}
PIPELINE_FORWARD_CONTEXT = dict(
    ignore_unknown_options=True,
    allow_extra_args=True,
    allow_interspersed_args=False,
)

pipeline_app = typer.Typer(
    help="Run and configure CapCruncher Snakemake workflows.",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
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


def should_touch_pipeline_outputs(options) -> bool:
    """Return whether a successful pipeline run should be followed by --touch."""
    return not (
        has_snakemake_option(options, "--dry-run", "-n")
        or has_snakemake_option(options, "--dryrun")
        or has_snakemake_option(options, "--touch")
    )


def get_capcruncher_config_dir() -> pathlib.Path:
    xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
    if xdg_config_home:
        return pathlib.Path(xdg_config_home).expanduser()

    return pathlib.Path.home() / ".config"


def get_pipeline_preset_dir() -> pathlib.Path:
    return get_capcruncher_config_dir() / "snakemake"


def normalise_pipeline_preset_name(preset: str) -> str:
    return LEGACY_PIPELINE_PRESET_ALIASES.get(preset, preset)


def resolve_pipeline_preset(preset: str) -> pathlib.Path:
    preset_path = pathlib.Path(preset).expanduser()
    if preset_path.exists():
        return preset_path.resolve()

    preset_name = normalise_pipeline_preset_name(preset)
    bundled_path = get_pipeline_preset_dir() / preset_name
    if bundled_path.exists():
        return bundled_path.resolve()

    raise click.ClickException(
        f"Unknown pipeline preset '{preset}'. Run 'capcruncher pipeline-init' to install presets or pass a profile path."
    )


def install_pipeline_preset(
    preset_name: str, output_dir: pathlib.Path, force: bool
) -> pathlib.Path:
    preset_name = normalise_pipeline_preset_name(preset_name)
    source_name = PIPELINE_PRESET_SOURCES[preset_name]
    source_dir = resources.files("capcruncher").joinpath(
        "pipeline", "profiles", source_name
    )
    destination_dir = output_dir / preset_name

    if destination_dir.exists() and not force:
        raise click.ClickException(
            f"Preset '{preset_name}' already exists at {destination_dir}. Use --force to overwrite it."
        )

    if destination_dir.exists() and force:
        shutil.rmtree(destination_dir)

    with resources.as_file(source_dir) as source_path:
        shutil.copytree(source_path, destination_dir)

    return destination_dir


def run_pipeline(
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
    elif should_touch_pipeline_outputs(pipeline_options):
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


def install_pipeline_presets(output_dir=None, preset_names=(), force=False):
    destination_root = output_dir or get_pipeline_preset_dir()
    destination_root.mkdir(parents=True, exist_ok=True)

    presets_to_install = preset_names or BUILTIN_PIPELINE_PRESETS
    installed = []
    for preset_name in presets_to_install:
        installed.append(install_pipeline_preset(preset_name, destination_root, force))

    click.echo(f"Installed {len(installed)} pipeline preset(s) to {destination_root}")
    for installed_preset in installed:
        click.echo(f"- {installed_preset.name}: {installed_preset}")


def _load_cookiecutter():
    try:
        from cookiecutter.main import cookiecutter
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.startswith("cookiecutter"):
            raise ModuleNotFoundError(
                "Cookiecutter is required to generate pipeline configuration. "
                "Install CapCruncher with the 'config' extra."
            ) from exc
        raise

    return cookiecutter


def configure_pipeline():
    import pathlib

    cookiecutter = _load_cookiecutter()
    fn = pathlib.Path(__file__).resolve()
    dir_cli = fn.parent
    dir_package = dir_cli.parent

    cookiecutter(str(dir_package / "pipeline" / "config"))


@pipeline_app.command("run", context_settings=PIPELINE_FORWARD_CONTEXT)
def pipeline_run(
    ctx: typer.Context,
    logo: bool = typer.Option(
        True,
        "--logo/--no-logo",
        help="Show the capcruncher logo.",
        show_default=True,
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        help="CapCruncher-managed execution preset name or a profile directory path.",
    ),
    scale_resources: float | None = typer.Option(
        None,
        "--scale-resources",
        help="Scale workflow memory and runtime requests.",
    ),
):
    """Run the CapCruncher Snakemake pipeline."""

    run_pipeline(
        tuple(ctx.args),
        logo=logo,
        preset=preset,
        scale_resources=scale_resources,
    )


@pipeline_app.command("init")
def pipeline_init_typer(
    output_dir: pathlib.Path | None = typer.Option(
        None,
        "--output-dir",
        file_okay=False,
        dir_okay=True,
        help="Directory where CapCruncher-managed pipeline presets should be installed.",
    ),
    preset_names: list[str] = typer.Option(
        None,
        "--preset",
        help="Install only the selected preset. Repeat to install multiple presets.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing preset directories if they already exist.",
    ),
):
    """Install CapCruncher-managed Snakemake presets."""

    invalid_presets = [
        preset_name
        for preset_name in (preset_names or [])
        if preset_name not in PIPELINE_PRESET_CHOICES
    ]
    if invalid_presets:
        raise click.ClickException(
            f"Unknown pipeline preset(s): {', '.join(invalid_presets)}"
        )

    install_pipeline_presets(output_dir, tuple(preset_names or ()), force)


@pipeline_app.command("config")
def pipeline_config_typer():
    """Create a CapCruncher pipeline configuration directory."""

    configure_pipeline()


def dispatch_pipeline_subcommand(pipeline_options) -> bool:
    if not pipeline_options or pipeline_options[0] not in PIPELINE_SUBCOMMANDS:
        return False

    pipeline_command = typer.main.get_command(pipeline_app)
    pipeline_command.main(
        args=list(pipeline_options),
        prog_name="capcruncher pipeline",
        standalone_mode=False,
    )
    return True


@cli.command(context_settings=PIPELINE_FORWARD_CONTEXT, name="pipeline")
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

    if dispatch_pipeline_subcommand(pipeline_options):
        return

    run_pipeline(
        pipeline_options,
        show_help=show_help,
        logo=logo,
        preset=preset,
        scale_resources=scale_resources,
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
    type=click.Choice(PIPELINE_PRESET_CHOICES),
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

    install_pipeline_presets(output_dir, preset_names, force)


@cli.command(name="pipeline-config", context_settings=CONTEXT_SETTINGS)
@click.version_option(get_capcruncher_version())
@click.option(
    "-i", "--input", "input_files", type=click.Path(exists=True), multiple=True
)
@click.option("--generate-design", is_flag=True)
def pipeline_config(*args, **kwargs):
    """Configures the data processing pipeline"""

    configure_pipeline()
