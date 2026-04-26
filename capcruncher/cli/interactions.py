from importlib import import_module
import os
import pathlib
import tempfile
from typing import Any, Literal

import typer

type FilePath = str | os.PathLike[str]


interactions_app = typer.Typer(
    help="Contains methods for interaction counting, storing, bedgraph generation, comparisons.",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)
compare_app = typer.Typer(
    help="Compare bedgraphs and CapCruncher cooler files.",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)


def _run_command(import_path: str, *args: Any, **kwargs: Any) -> None:
    module_name, function_name = import_path.rsplit(":", 1)
    module = import_module(module_name)
    getattr(module, function_name)(*args, **kwargs)


def _valid_viewpoint_names(viewpoint_path: FilePath) -> list[str]:
    import pandas as pd

    viewpoints = pd.read_csv(
        viewpoint_path,
        sep="\t",
        header=None,
        usecols=[3],
        names=["name"],
    )
    return viewpoints["name"].dropna().astype(str).drop_duplicates().tolist()


def _parquet_files(path: FilePath) -> list[pathlib.Path]:
    path = pathlib.Path(path)
    if path.is_dir():
        return sorted(path.glob("*.parquet"))
    return [path]


def _write_countable_reporters(
    reporters: FilePath, viewpoint_path: FilePath, output_dir: FilePath
) -> pathlib.Path:
    import pandas as pd

    valid_viewpoints = _valid_viewpoint_names(viewpoint_path)
    if not valid_viewpoints:
        raise ValueError(f"No viewpoints found in {viewpoint_path}")

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for index, parquet_file in enumerate(_parquet_files(reporters)):
        reporters_df = pd.read_parquet(parquet_file)
        reporters_df["viewpoint"] = reporters_df["viewpoint"].astype(str)
        invalid_viewpoints = sorted(
            set(reporters_df["viewpoint"].dropna()) - set(valid_viewpoints)
        )
        if invalid_viewpoints:
            raise ValueError(
                "Reporter file contains viewpoint values not present in "
                f"{viewpoint_path}: {invalid_viewpoints}"
            )

        for column in reporters_df.select_dtypes(include="category").columns:
            reporters_df[column] = reporters_df[column].cat.remove_unused_categories()

        reporters_df["viewpoint"] = pd.Categorical(
            reporters_df["viewpoint"],
            categories=valid_viewpoints,
        )
        reporters_df.to_parquet(output_dir / f"part-{index}.parquet", index=False)

    return output_dir


def count_interactions(
    reporters: FilePath,
    output: FilePath = "CC_cooler.hdf5",
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    fragment_map: FilePath | None = None,
    viewpoint_path: FilePath | None = None,
    n_cores: int = 1,
    assay: Literal["capture", "tri", "tiled"] = "capture",
    executor: Literal["local", "process", "ray"] = "local",
    **kwargs: Any,
) -> FilePath:
    from capcruncher_tools.api import count_interactions as count_interactions_records

    with tempfile.TemporaryDirectory() as tmpdir:
        countable_reporters = _write_countable_reporters(
            reporters=reporters,
            viewpoint_path=viewpoint_path,
            output_dir=tmpdir,
        )

        clr = count_interactions_records(
            reporters=countable_reporters,
            output=output,
            remove_exclusions=remove_exclusions,
            remove_viewpoint=remove_viewpoint,
            subsample=subsample,
            fragment_map=fragment_map,
            viewpoint_path=viewpoint_path,
            n_cores=n_cores,
            assay=assay,
            executor=executor,
            **kwargs,
        )

    return clr


@interactions_app.callback()
def interactions() -> None:
    """Contains methods for interaction counting, storing, bedgraph generation, comparisons."""


@interactions_app.command()
def deduplicate(
    slices: str = typer.Argument(...),
    output: str = typer.Option(
        "deduplicated_slices/",
        "-o",
        "--output",
        help="Output prefix for directory of deduplicated slices.",
    ),
    statistics: str = typer.Option(
        "",
        "--statistics",
        help="Output prefix for stats file(s).",
    ),
    sample_name: str = typer.Option(
        "sample",
        "--sample-name",
        help="Name of sample e.g. DOX_treated_1.",
    ),
    read_type: Literal["flashed", "pe"] = typer.Option(
        "flashed",
        "--read-type",
        help="Type of read.",
    ),
) -> None:
    """
    Identifies and removes duplicated aligned fragments.
    """
    _run_command(
        "capcruncher.api.interactions_deduplicate:deduplicate",
        slices=slices,
        output=output,
        statistics=statistics,
        sample_name=sample_name,
        read_type=read_type,
    )


@interactions_app.command()
def pileup(
    uri: str = typer.Argument(...),
    viewpoint_names: list[str] | None = typer.Option(
        None,
        "-n",
        "--viewpoint-names",
        "--viewpoint_names",
        help="Viewpoint to extract and convert to bedgraph. If not provided, transform all.",
    ),
    output_prefix: str = typer.Option(
        "",
        "-o",
        "--output-prefix",
        "--output_prefix",
        help="Output prefix for bedgraphs.",
    ),
    normalisation: Literal["raw", "n_cis", "region"] = typer.Option(
        "raw",
        "--normalisation",
        help="Method to use interaction normalisation.",
    ),
    normalisation_regions: str | None = typer.Option(
        None,
        "--normalisation-regions",
        help="Regions to use for interaction normalisation. The --normalisation method MUST be 'region'.",
    ),
    binsize: int = typer.Option(
        0,
        "--binsize",
        help="Binsize to use for converting bedgraph to evenly sized genomic bins.",
    ),
    gzip: bool = typer.Option(
        False,
        "--gzip",
        help="Compress output using gzip.",
    ),
    scale_factor: float = typer.Option(
        1e6,
        "--scale-factor",
        help="Scale factor to use for bedgraph normalisation.",
    ),
    sparse: bool = typer.Option(
        True,
        "--sparse/--dense",
        help="Produce bedgraph containing just positive bins (sparse) or all bins (dense).",
    ),
    format: Literal["bedgraph", "bigwig"] = typer.Option(
        "bedgraph",
        "-f",
        "--format",
        help="Output file format.",
    ),
) -> None:
    """
    Extract reporters from a capture experiment and generate a bedgraph or bigWig file.
    """
    _run_command(
        "capcruncher.api.pileup:pileup",
        uri=uri,
        viewpoint_names=viewpoint_names,
        output_prefix=output_prefix,
        normalisation=normalisation,
        normalisation_regions=normalisation_regions,
        binsize=binsize,
        gzip=gzip,
        scale_factor=scale_factor,
        sparse=sparse,
        format=format,
    )


@interactions_app.command()
def count(
    reporters: str = typer.Argument(...),
    output: str = typer.Option(
        "CC_cooler.hdf5",
        "-o",
        "--output",
        help="Name of output file.",
    ),
    remove_exclusions: bool = typer.Option(
        False,
        "--remove-exclusions",
        "--remove_exclusions",
        help="Prevents analysis of fragments marked as proximity exclusions.",
    ),
    remove_capture: bool = typer.Option(
        False,
        "--remove-capture",
        "--remove_capture",
        help="Prevents analysis of capture fragment interactions.",
    ),
    subsample: float = typer.Option(
        0,
        "--subsample",
        help="Subsamples reporters before analysis of interactions.",
    ),
    fragment_map: str | None = typer.Option(
        None,
        "-f",
        "--fragment-map",
        help="Path to digested genome bed file.",
    ),
    viewpoint_path: str | None = typer.Option(
        None,
        "-v",
        "--viewpoint-path",
        help="Path to viewpoints file.",
    ),
    n_cores: int = typer.Option(
        1,
        "-p",
        "--n-cores",
        "--n_cores",
        help="Number of cores to use for counting.",
    ),
    assay: Literal["capture", "tri", "tiled"] = typer.Option(
        "capture",
        "--assay",
    ),
    executor: Literal["local", "process", "ray"] = typer.Option(
        "local",
        "--executor",
        help="Runtime used for per-viewpoint counting.",
    ),
) -> None:
    """
    Determines the number of captured restriction fragment interactions genome wide.
    """
    count_interactions(
        reporters=reporters,
        output=output,
        remove_exclusions=remove_exclusions,
        remove_viewpoint=remove_capture,
        subsample=subsample,
        fragment_map=fragment_map,
        viewpoint_path=viewpoint_path,
        n_cores=n_cores,
        assay=assay,
        executor=executor,
    )


def store_fragments(
    counts: str = typer.Argument(...),
    fragment_map: str = typer.Option(
        ...,
        "-f",
        "--fragment-map",
        help="Path to digested genome bed file.",
    ),
    viewpoint_path: str = typer.Option(
        ...,
        "-v",
        "--viewpoint-path",
        help="Path to viewpoints file.",
    ),
    viewpoint_name: str = typer.Option(
        "",
        "-n",
        "--viewpoint-name",
        help="Name of viewpoint to store.",
    ),
    genome: str = typer.Option(
        "",
        "-g",
        "--genome",
        help="Name of genome.",
    ),
    suffix: str = typer.Option(
        "",
        "--suffix",
        help="Suffix to append after the capture name for the output file.",
    ),
    output: str = typer.Option(
        "out.hdf5",
        "-o",
        "--output",
        help="Name of output file. (Cooler formatted hdf5 file).",
    ),
) -> None:
    """
    Stores restriction fragment interaction combinations at the restriction fragment level.
    """
    _run_command(
        "capcruncher.api.storage:fragments",
        counts=counts,
        fragment_map=fragment_map,
        output=output,
        viewpoint_path=viewpoint_path,
        viewpoint_name=viewpoint_name,
        genome=genome,
        suffix=suffix,
    )


def store_bins(
    cooler_path: str = typer.Argument(...),
    binsizes: list[int] = typer.Option(
        [5000],
        "-b",
        "--binsizes",
        help="Binsizes to use for windowing.",
    ),
    normalise: bool = typer.Option(
        False,
        "--normalise",
        help="Enables normalisation of interaction counts during windowing.",
    ),
    overlap_fraction: float = typer.Option(
        0.5,
        "--overlap-fraction",
        "--overlap_fraction",
        help="Minimum overlap between genomic bins and restriction fragments for overlap.",
    ),
    n_cores: int = typer.Option(
        4,
        "-p",
        "--n-cores",
        "--n_cores",
        help="Number of cores used for binning.",
    ),
    scale_factor: float = typer.Option(
        1e6,
        "--scale-factor",
        help="Scaling factor used for normalisation.",
    ),
    conversion_tables: str | None = typer.Option(
        None,
        "--conversion-tables",
        "--conversion_tables",
        help="Pickle file containing pre-computed fragment -> bin conversions.",
    ),
    output: str = typer.Option(
        "out.hdf5",
        "-o",
        "--output",
        help="Name of output file. (Cooler formatted hdf5 file).",
    ),
    assay: Literal["capture", "tri", "tiled"] = typer.Option(
        "capture",
        "--assay",
    ),
) -> None:
    """
    Convert a cooler group containing restriction fragments to constant genomic windows.
    """
    _run_command(
        "capcruncher.api.storage:bins",
        cooler_path=cooler_path,
        output=output,
        binsizes=tuple(binsizes),
        normalise=normalise,
        scale_factor=scale_factor,
        overlap_fraction=overlap_fraction,
        conversion_tables=conversion_tables,
        n_cores=n_cores,
        assay=assay,
    )


@interactions_app.command(name="merge")
def store_merge(
    coolers: list[str] = typer.Argument(...),
    output: str = typer.Option(
        ...,
        "-o",
        "--output",
        help="Output file name.",
    ),
) -> None:
    """
    Merges CapCruncher HDF5 files together.
    """
    _run_command(
        "capcruncher.api.storage:merge_coolers",
        coolers=tuple(coolers),
        output=output,
    )


@compare_app.callback()
def compare() -> None:
    """Compare bedgraphs and CapCruncher cooler files."""


@compare_app.command(name="concat")
def bedgraphs_concat(
    infiles: list[str] = typer.Argument(...),
    format: Literal["auto", "bedgraph", "cooler"] = typer.Option(
        "cooler",
        "-f",
        "--format",
        help="Input file format.",
    ),
    output: str = typer.Option(
        "union.tsv",
        "-o",
        "--output",
        help="Output file name.",
    ),
    viewpoint: str | None = typer.Option(
        None,
        "-v",
        "--viewpoint",
        help="Viewpoint to extract.",
    ),
    resolution: int | None = typer.Option(
        None,
        "-r",
        "--resolution",
        help="Resolution to extract.",
    ),
    region: str | None = typer.Option(
        None,
        "--region",
        help="Limit to specific coordinates in the format chrom:start-end.",
    ),
    normalisation: Literal["raw", "n_cis", "region"] = typer.Option(
        "raw",
        "--normalisation",
        help="Method to use interaction normalisation.",
    ),
    normalisation_regions: str | None = typer.Option(
        None,
        "--normalisation-regions",
        help="Regions to use for interaction normalisation. The --normalisation method MUST be 'region'.",
    ),
    scale_factor: float = typer.Option(
        1e6,
        "--scale-factor",
        "--scale_factor",
        help="Scale factor to use for bedgraph normalisation.",
    ),
    n_cores: int = typer.Option(
        1,
        "-p",
        "--n-cores",
        "--n_cores",
        help="Number of cores to use for extracting bedgraphs.",
    ),
) -> None:
    _run_command(
        "capcruncher.api.interactions_compare:concat",
        infiles=tuple(infiles),
        format=format,
        output=output,
        viewpoint=viewpoint,
        resolution=resolution,
        region=region,
        normalisation=normalisation,
        normalisation_regions=normalisation_regions,
        scale_factor=scale_factor,
        n_cores=n_cores,
    )


@compare_app.command(name="summarise")
def bedgraphs_summarise(
    infile: str = typer.Argument(...),
    design_matrix: str | None = typer.Option(
        None,
        "-d",
        "--design-matrix",
        help="Design matrix file.",
    ),
    output_prefix: str | None = typer.Option(
        None,
        "-o",
        "--output-prefix",
        help="Output file prefix.",
    ),
    output_format: Literal["bedgraph", "tsv"] = typer.Option(
        "bedgraph",
        "-f",
        "--output-format",
        help="Output file format.",
    ),
    summary_methods: list[str] | None = typer.Option(
        None,
        "-m",
        "--summary-methods",
        help="Summary methods to use for aggregation.",
    ),
    group_names: list[str] | None = typer.Option(
        None,
        "-n",
        "--group-names",
        help="Group names for aggregation.",
    ),
    group_columns: list[str] | None = typer.Option(
        None,
        "-c",
        "--group-columns",
        help="Column names/numbers for aggregation.",
    ),
    perform_subtractions: bool = typer.Option(
        False,
        "--subtraction",
        help="Perform subtraction between aggregated groups.",
    ),
    suffix: str = typer.Option(
        "",
        "--suffix",
        help="Add a suffix before the file extension.",
    ),
) -> None:
    _run_command(
        "capcruncher.api.interactions_compare:summarise",
        infile=infile,
        design_matrix=design_matrix,
        output_prefix=output_prefix,
        output_format=output_format,
        summary_methods=tuple(summary_methods or ()),
        group_names=tuple(group_names or ()),
        group_columns=tuple(group_columns or ()),
        suffix=suffix,
        perform_subtractions=perform_subtractions,
    )


def _run_differential(
    interaction_files: list[str],
    output_prefix: str,
    viewpoint: str,
    design_matrix: str,
    contrast: str,
    regions_of_interest: str | None,
    viewpoint_distance: int | None,
    threshold_count: float,
    threshold_q: float,
) -> None:
    _run_command(
        "capcruncher.api.interactions_differential:differential",
        interaction_files=tuple(interaction_files),
        output_prefix=output_prefix,
        viewpoint=viewpoint,
        design_matrix=design_matrix,
        contrast=contrast,
        regions_of_interest=regions_of_interest,
        viewpoint_distance=viewpoint_distance,
        threshold_count=threshold_count,
        threshold_q=threshold_q,
    )


def bedgraphs_differential(
    interaction_files: list[str] = typer.Argument(...),
    output_prefix: str = typer.Option(
        "differential",
        "-o",
        "--output-prefix",
        help="Output file prefix.",
    ),
    viewpoint: str = typer.Option(
        ...,
        "-v",
        "--viewpoint",
        help="Viewpoint to extract.",
    ),
    design_matrix: str = typer.Option(
        ...,
        "-d",
        "--design-matrix",
        help="Design matrix file.",
    ),
    contrast: str = typer.Option(
        "condition",
        "-c",
        "--contrast",
        help="Contrast to test.",
    ),
    regions_of_interest: str | None = typer.Option(
        None,
        "-r",
        "--regions-of-interest",
        help="Regions of interest to test for differential interactions.",
    ),
    viewpoint_distance: int | None = typer.Option(
        None,
        "--viewpoint-distance",
        help="Distance from viewpoint to test for differential interactions.",
    ),
    threshold_count: float = typer.Option(
        20,
        "--threshold-count",
        help="Minimum number of interactions to test for differential interactions.",
    ),
    threshold_q: float = typer.Option(
        0.05,
        "--threshold-q",
        help="Minimum q-value to test for differential interactions.",
    ),
) -> None:
    """Perform differential testing on CapCruncher HDF5 files."""
    _run_differential(
        interaction_files=interaction_files,
        output_prefix=output_prefix,
        viewpoint=viewpoint,
        design_matrix=design_matrix,
        contrast=contrast,
        regions_of_interest=regions_of_interest,
        viewpoint_distance=viewpoint_distance,
        threshold_count=threshold_count,
        threshold_q=threshold_q,
    )


interactions_app.command(name="counts-to-cooler")(store_fragments)
interactions_app.command(name="fragments-to-bins")(store_bins)
interactions_app.command(name="bin")(store_bins)
interactions_app.command(name="differential")(bedgraphs_differential)
compare_app.command(name="differential")(bedgraphs_differential)
interactions_app.add_typer(compare_app, name="compare")

cli = typer.main.get_command(interactions_app)
