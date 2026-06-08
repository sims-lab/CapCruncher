import os
import pathlib
import shutil
import subprocess
from collections.abc import Sequence
from importlib import resources

import typer

from capcruncher.cli.common import HELP_SETTINGS
from capcruncher.dependencies import DependencyVersionError, require_capcruncher_tools

type PipelineOptions = Sequence[str]


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
PIPELINE_FORWARD_CONTEXT = dict(
    ignore_unknown_options=True,
    allow_extra_args=True,
    allow_interspersed_args=False,
)

pipeline_app = typer.Typer(
    help="Run and configure CapCruncher Snakemake workflows.",
    context_settings=HELP_SETTINGS,
    no_args_is_help=True,
)


def has_snakemake_option(
    options: PipelineOptions, long_name: str, short_name: str | None = None
) -> bool:
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

    raise typer.BadParameter(
        f"Unknown pipeline preset '{preset}'. Run 'capcruncher pipeline init' to install presets or pass a profile path."
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
        raise typer.BadParameter(
            f"Preset '{preset_name}' already exists at {destination_dir}. Use --force to overwrite it."
        )

    if destination_dir.exists() and force:
        shutil.rmtree(destination_dir)

    with resources.as_file(source_dir) as source_path:
        shutil.copytree(source_path, destination_dir)

    return destination_dir


def run_pipeline(
    pipeline_options: PipelineOptions,
    show_help: bool = False,
    logo: bool = True,
    preset: str | None = None,
    scale_resources: float | None = None,
) -> None:
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
        cmd.append("--help")
        _completed = subprocess.run(cmd, capture_output=True, shell=False, text=True)
        _print_pipeline_run_help(_completed.stdout)
        raise typer.Exit()

    if pipeline_options:
        excluded_options = ["--version", "make", "run", "show"]

        cmd.extend(
            [option for option in pipeline_options if option not in excluded_options]
        )

    if preset:
        if has_snakemake_option(pipeline_options, "--profile"):
            typer.echo("Use either --preset or --profile, not both.", err=True)
            raise typer.Exit(2)
        cmd.extend(["--profile", str(resolve_pipeline_preset(preset))])

    try:
        require_capcruncher_tools()
    except DependencyVersionError as exc:
        typer.secho(str(exc), err=True, fg=typer.colors.RED)
        raise typer.Exit(1) from exc

    # Implicitly deal with a missing --cores option
    if not has_snakemake_option(pipeline_options, "--cores", "-c"):
        cmd.extend(["--cores", "1"])

    # Add the --show-failed-logs option if it is not already present
    if not has_snakemake_option(pipeline_options, "--show-failed-logs"):
        cmd.append("--show-failed-logs")

    if logo:
        with open(dir_package / "data" / "logo.txt", encoding="utf-8") as f:
            typer.echo(f.read())

    env = os.environ.copy()
    if scale_resources is not None:
        env["SCALE_RESOURCES"] = str(scale_resources)

    # Run the pipeline
    _completed = subprocess.run(cmd, env=env)

    # If the pipeline fails, exit with the return code
    if _completed.returncode != 0:
        raise typer.Exit(_completed.returncode)

    if should_touch_pipeline_outputs(pipeline_options):
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


def install_pipeline_presets(
    output_dir: pathlib.Path | None = None,
    preset_names: Sequence[str] = (),
    force: bool = False,
) -> None:
    destination_root = output_dir or get_pipeline_preset_dir()
    destination_root.mkdir(parents=True, exist_ok=True)

    presets_to_install = preset_names or BUILTIN_PIPELINE_PRESETS
    installed = []
    for preset_name in presets_to_install:
        installed.append(install_pipeline_preset(preset_name, destination_root, force))

    typer.echo(f"Installed {len(installed)} pipeline preset(s) to {destination_root}")
    for installed_preset in installed:
        typer.echo(f"- {installed_preset.name}: {installed_preset}")


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
    cookiecutter = _load_cookiecutter()
    fn = pathlib.Path(__file__).resolve()
    dir_cli = fn.parent
    dir_package = dir_cli.parent

    cookiecutter(str(dir_package / "pipeline" / "config"))


def _option_value(options: PipelineOptions, index: int, option: str) -> tuple[str, int]:
    if index + 1 >= len(options):
        raise typer.BadParameter(f"Option '{option}' requires a value.")
    return options[index + 1], index + 1


def _float_option_value(value: str, option: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise typer.BadParameter(
            f"Option '{option}' requires a numeric value."
        ) from exc


def _print_pipeline_run_help(snakemake_help: str) -> None:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table

    console = Console()
    console.print(
        Panel.fit(
            "\n".join(
                [
                    "[bold]capcruncher pipeline run[/bold] "
                    "[CAPCRUNCHER OPTIONS] [SNAKEMAKE OPTIONS] [TARGETS]",
                    "",
                    "CapCruncher consumes the options listed below. Everything else is passed to Snakemake.",
                ]
            ),
            title="Usage",
            border_style="cyan",
        )
    )

    capcruncher_options = Table(
        title="CapCruncher Options",
        show_header=True,
        header_style="bold cyan",
    )
    capcruncher_options.add_column("Option", no_wrap=True)
    capcruncher_options.add_column("Description")
    capcruncher_options.add_row(
        "--preset TEXT",
        "Use an installed CapCruncher Snakemake profile. "
        "Aliases such as 'local' are accepted. Cannot be combined with --profile.",
    )
    capcruncher_options.add_row(
        "--scale-resources FLOAT",
        "Set CapCruncher's SCALE_RESOURCES environment value for workflow resource functions.",
    )
    capcruncher_options.add_row(
        "--logo / --no-logo",
        "Show or suppress the CapCruncher logo before running Snakemake.",
    )
    capcruncher_options.add_row(
        "-h, --help",
        "Show this combined CapCruncher and Snakemake help.",
    )
    console.print(capcruncher_options)

    console.print(
        Panel.fit(
            "\n".join(
                [
                    "[bold]Snakemake options and targets[/bold]",
                    "",
                    "Use Snakemake options after the CapCruncher options, for example:",
                    "[green]capcruncher pipeline run --preset local -c 4 -n results/file.parquet[/green]",
                    "",
                    "Common Snakemake options include -c/--cores, -n/--dry-run, "
                    "--config, --configfile, --profile, and explicit workflow targets.",
                ]
            ),
            border_style="green",
        )
    )

    output = snakemake_help.replace("usage: snakemake", "Snakemake usage: snakemake")
    console.print("[bold]Underlying Snakemake Help[/bold]")
    console.print(output.rstrip())


def _parse_pipeline_run_options(
    options: PipelineOptions,
) -> tuple[tuple[str, ...], bool, bool, str | None, float | None]:
    remaining: list[str] = []
    logo = True
    preset = None
    scale_resources = None
    show_help = False
    index = 0

    while index < len(options):
        option = options[index]

        if option in {"-h", "--help"}:
            show_help = True
        elif option == "--logo":
            logo = True
        elif option == "--no-logo":
            logo = False
        elif option == "--preset":
            preset, index = _option_value(options, index, option)
        elif option.startswith("--preset="):
            preset = option.split("=", 1)[1]
        elif option == "--scale-resources":
            value, index = _option_value(options, index, option)
            scale_resources = _float_option_value(value, option)
        elif option.startswith("--scale-resources="):
            scale_resources = _float_option_value(
                option.split("=", 1)[1], "--scale-resources"
            )
        else:
            remaining.append(option)

        index += 1

    return tuple(remaining), show_help, logo, preset, scale_resources


def _run_pipeline_init(
    output_dir: pathlib.Path | None,
    preset_names: list[str] | None,
    force: bool,
) -> None:
    """Install CapCruncher-managed Snakemake presets."""

    invalid_presets = [
        preset_name
        for preset_name in (preset_names or [])
        if preset_name not in PIPELINE_PRESET_CHOICES
    ]
    if invalid_presets:
        raise typer.BadParameter(
            f"Unknown pipeline preset(s): {', '.join(invalid_presets)}"
        )

    install_pipeline_presets(output_dir, tuple(preset_names or ()), force)


@pipeline_app.callback()
def pipeline(ctx: typer.Context) -> None:
    """Run and configure CapCruncher Snakemake workflows."""


@pipeline_app.command(
    name="run",
    context_settings={
        **PIPELINE_FORWARD_CONTEXT,
        "help_option_names": [],
    },
)
def pipeline_run(ctx: typer.Context) -> None:
    """Run the data processing pipeline."""

    run_options, show_help, logo, preset, scale_resources = _parse_pipeline_run_options(
        tuple(ctx.args)
    )
    run_pipeline(
        run_options,
        show_help=show_help,
        logo=logo,
        preset=preset,
        scale_resources=scale_resources,
    )


@pipeline_app.command(name="init")
def pipeline_init_command(
    output_dir: pathlib.Path | None = typer.Option(
        None,
        "--output-dir",
        file_okay=False,
        dir_okay=True,
        help="Directory where CapCruncher-managed pipeline presets should be installed.",
    ),
    preset_names: list[str] | None = typer.Option(
        None,
        "--preset",
        help="Install only the selected preset. Repeat to install multiple presets.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing preset directories if they already exist.",
    ),
) -> None:
    """Install CapCruncher-managed Snakemake presets."""

    _run_pipeline_init(output_dir, preset_names, force)


@pipeline_app.command(name="config")
def pipeline_config_command(
    list_profiles: bool = typer.Option(
        False,
        "--list-profiles",
        help="Print stored genome profiles and exit.",
    ),
) -> None:
    """Configure the data processing pipeline."""
    if list_profiles:
        from capcruncher.cli.genome import genome_profile_list

        genome_profile_list()
        raise typer.Exit()

    configure_pipeline()


@pipeline_app.command(name="design")
def pipeline_design_command(
    output: pathlib.Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Write inferred design matrix to this TSV path instead of printing.",
    ),
    condition_pattern: str | None = typer.Option(
        None,
        "--condition-pattern",
        help=(
            "Regex with a named group 'condition' for custom filename conventions. "
            "Default: everything before the last underscore is the condition."
        ),
    ),
) -> None:
    """Infer a design matrix from *.fastq.gz files in the current directory.

    Expected filename convention: <CONDITION>_<REPLICATE>_R[12].fastq[.gz]

    Prints a preview table; use --output/-o to save as TSV.
    """
    import glob

    from capcruncher.pipeline.utils import infer_design_from_fastqs

    fastqs = sorted(glob.glob("*.fastq.gz") + glob.glob("*.fastq"))
    if not fastqs:
        typer.secho(
            "No *.fastq[.gz] files found in current directory.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(1)

    design = infer_design_from_fastqs(fastqs, condition_pattern=condition_pattern)

    try:
        from rich.console import Console
        from rich.table import Table

        table = Table(title="Inferred Design Matrix", show_lines=True)
        for col in design.columns:
            table.add_column(col, style="cyan")
        for row in design.itertuples(index=False):
            table.add_row(*[str(v) if v is not None else "" for v in row])
        Console().print(table)
    except ImportError:
        typer.echo(design.to_string(index=False))

    if output:
        design.to_csv(output, sep="\t", index=False)
        typer.echo(f"Design matrix written to {output}")


def pipeline_init(
    output_dir: pathlib.Path | None = typer.Option(
        None,
        "--output-dir",
        file_okay=False,
        dir_okay=True,
        help="Directory where CapCruncher-managed pipeline presets should be installed.",
    ),
    preset_names: list[str] | None = typer.Option(
        None,
        "--preset",
        help="Install only the selected preset. Repeat to install multiple presets.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Overwrite existing preset directories if they already exist.",
    ),
) -> None:
    """Installs CapCruncher-managed Snakemake presets."""

    typer.echo(
        "Warning: 'capcruncher pipeline-init' is deprecated. "
        "Use 'capcruncher pipeline init ...' instead.",
        err=True,
    )

    invalid_presets = [
        preset_name
        for preset_name in (preset_names or [])
        if preset_name not in PIPELINE_PRESET_CHOICES
    ]
    if invalid_presets:
        raise typer.BadParameter(
            f"Unknown pipeline preset(s): {', '.join(invalid_presets)}"
        )

    _run_pipeline_init(output_dir, preset_names, force)


def pipeline_config(
    input_files: list[pathlib.Path] | None = typer.Option(
        None,
        "-i",
        "--input",
        exists=True,
        help="Input files.",
    ),
    generate_design: bool = typer.Option(
        False,
        "--generate-design",
        help="Generate a design matrix.",
    ),
) -> None:
    """Configures the data processing pipeline"""

    typer.echo(
        "Warning: 'capcruncher pipeline-config' is deprecated. "
        "Use 'capcruncher pipeline config ...' instead.",
        err=True,
    )
    configure_pipeline()


cli = typer.main.get_command(pipeline_app)
