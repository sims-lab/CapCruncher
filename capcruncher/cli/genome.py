import os
import pathlib

import typer
import yaml

from capcruncher.cli.common import HELP_SETTINGS

genome_app = typer.Typer(
    help="Contains methods for genome digestion and genome profile management.",
    context_settings=HELP_SETTINGS,
    no_args_is_help=True,
)


def _genome_profiles_dir() -> pathlib.Path:
    xdg = os.environ.get("XDG_CONFIG_HOME")
    base = (
        pathlib.Path(xdg).expanduser() if xdg else pathlib.Path.home() / ".capcruncher"
    )
    return base / "genomes"


@genome_app.callback()
def genome():
    """Contains methods for genome digestion."""


@genome_app.command()
def digest(
    input_fasta: str = typer.Argument(...),
    recognition_site: str = typer.Option(
        ...,
        "-r",
        "--recognition-site",
        "--recognition_site",
        help="Recognition enzyme or sequence.",
    ),
    logfile: str = typer.Option(
        "genome_digest.log",
        "-l",
        "--logfile",
        help="Path for digestion log file.",
    ),
    output_file: str = typer.Option(
        "genome_digested.bed",
        "-o",
        "--output-file",
        "--output_file",
        help="Output file path.",
    ),
    remove_cutsite: bool = typer.Option(
        True,
        "--remove-cutsite/--keep-cutsite",
        "--remove_cutsite/--keep_cutsite",
        help="Exclude the recognition sequence from the output.",
    ),
    sort: bool = typer.Option(
        False,
        "--sort",
        help="Sorts the output bed file by chromosome and start coord.",
    ),
):
    """
    Performs in silico digestion of a genome in fasta format.

    Digests the supplied genome fasta file and generates a bed file containing the
    locations of all restriction fragments produced by the supplied restriction enzyme.

    A log file recording the number of restriction fragments for the suplied genome is also
    generated.
    """
    from capcruncher.api.genome import digest_genome

    digest_genome(
        input_fasta=input_fasta,
        recognition_site=recognition_site,
        output_file=output_file,
        logfile=logfile,
        remove_cutsite=remove_cutsite,
        sort=sort,
    )


@genome_app.command(name="add")
def genome_profile_add(
    name: str = typer.Argument(..., help="Profile name (e.g. hg38, mm10)."),
    fasta: str = typer.Option(..., prompt=True, help="Path to genome FASTA."),
    aligner_index: str = typer.Option(
        ..., prompt=True, help="Path to aligner index prefix."
    ),
    chrom_sizes: str = typer.Option(..., prompt=True, help="Path to chrom.sizes file."),
    organism: str = typer.Option(
        "", prompt=True, help="Organism name (e.g. 'Homo sapiens')."
    ),
    twobit: str = typer.Option(
        "", prompt=True, help="Path to .2bit file (optional, for custom hub genomes)."
    ),
    custom: bool = typer.Option(
        False, prompt=True, help="Custom genome (not in UCSC database)?"
    ),
) -> None:
    """Save a genome profile for reuse in pipeline configs."""
    profiles_dir = _genome_profiles_dir()
    profiles_dir.mkdir(parents=True, exist_ok=True)
    profile = {
        "name": name,
        "organism": organism or None,
        "fasta": fasta,
        "aligner_index": aligner_index,
        "chrom_sizes": chrom_sizes,
        "twobit": twobit or None,
        "custom": custom,
    }
    dest = profiles_dir / f"{name}.yml"
    dest.write_text(yaml.dump(profile, default_flow_style=False))
    typer.secho(f"Genome profile '{name}' saved to {dest}", fg=typer.colors.GREEN)


@genome_app.command(name="list")
def genome_profile_list() -> None:
    """List all stored genome profiles."""
    profiles_dir = _genome_profiles_dir()
    profiles = sorted(profiles_dir.glob("*.yml")) if profiles_dir.exists() else []
    if not profiles:
        typer.echo(
            "No genome profiles found. Use `capcruncher genome add` to create one."
        )
        return
    try:
        from rich.console import Console
        from rich.table import Table

        table = Table(title="Genome Profiles", show_lines=True)
        table.add_column("Name", style="cyan")
        table.add_column("Organism")
        table.add_column("FASTA")
        for p in profiles:
            data = yaml.safe_load(p.read_text())
            table.add_row(
                data.get("name", p.stem),
                data.get("organism") or "",
                data.get("fasta", ""),
            )
        Console().print(table)
    except ImportError:
        for p in profiles:
            data = yaml.safe_load(p.read_text())
            typer.echo(
                f"{data.get('name', p.stem)}\t{data.get('organism') or ''}\t{data.get('fasta', '')}"
            )


@genome_app.command(name="show")
def genome_profile_show(
    name: str = typer.Argument(..., help="Profile name to display."),
) -> None:
    """Print the full YAML for a stored genome profile."""
    profile_path = _genome_profiles_dir() / f"{name}.yml"
    if not profile_path.exists():
        typer.secho(f"Profile '{name}' not found.", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    typer.echo(profile_path.read_text())


@genome_app.command(name="remove")
def genome_profile_remove(
    name: str = typer.Argument(..., help="Profile name to delete."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation prompt."),
) -> None:
    """Delete a stored genome profile."""
    profile_path = _genome_profiles_dir() / f"{name}.yml"
    if not profile_path.exists():
        typer.secho(f"Profile '{name}' not found.", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    if not yes:
        typer.confirm(f"Delete profile '{name}'?", abort=True)
    profile_path.unlink()
    typer.secho(f"Profile '{name}' removed.", fg=typer.colors.YELLOW)


cli = typer.main.get_command(genome_app)
