from importlib import import_module
from typing import Annotated, Any

import typer

HELP_SETTINGS = {"help_option_names": ["-h", "--help"]}

CompressionLevelOption = Annotated[
    int,
    typer.Option(
        "--compression-level",
        "--compression_level",
        min=0,
        max=9,
        help="Level of compression for output files.",
    ),
]
MinimumSliceLengthOption = Annotated[
    int,
    typer.Option(
        "--minimum-slice-length",
        "--minimum_slice_length",
        min=1,
    ),
]
NCoresOption = Annotated[
    int,
    typer.Option(
        "-p",
        "--n-cores",
        "--n_cores",
        min=1,
        help="Number of cores to use.",
    ),
]
NReadsOption = Annotated[
    int,
    typer.Option(
        "-n",
        "--n-reads",
        "--n_reads",
        min=1,
        help="Number of reads per fastq file.",
    ),
]
SubsampleOption = Annotated[
    float,
    typer.Option(
        "--subsample",
        min=0,
        max=1,
        help="Subsamples reporters before analysis of interactions.",
    ),
]


def run_imported(import_path: str, *args: Any, **kwargs: Any) -> None:
    """Import and run a command implementation on demand."""
    module_name, function_name = import_path.rsplit(":", 1)
    module = import_module(module_name)
    getattr(module, function_name)(*args, **kwargs)
