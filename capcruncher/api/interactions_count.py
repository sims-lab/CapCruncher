import os
import pathlib
import tempfile
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, PositiveInt, field_validator

type FilePath = str | os.PathLike[str]


class InteractionCountOptions(BaseModel):
    """Validated options for reporter interaction counting."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    reporters: FilePath
    output: FilePath = "CC_cooler.hdf5"
    remove_exclusions: bool = False
    remove_viewpoint: bool = False
    subsample: float = Field(default=0, ge=0, le=1)
    fragment_map: FilePath | None = None
    viewpoint_path: FilePath
    n_cores: PositiveInt = 1
    assay: Literal["capture", "tri", "tiled"] = "capture"
    executor: Literal["local", "process", "ray"] = "local"

    @field_validator("reporters", "fragment_map", "viewpoint_path")
    @classmethod
    def existing_input_path(cls, value: FilePath | None) -> FilePath | None:
        if value is not None and not pathlib.Path(value).exists():
            raise ValueError(f"Input path does not exist: {value}")
        return value


def valid_viewpoint_names(viewpoint_path: FilePath) -> list[str]:
    """Return unique viewpoint names from a BED-like viewpoint file."""
    import pandas as pd

    viewpoints = pd.read_csv(
        viewpoint_path,
        sep="\t",
        header=None,
        usecols=[3],
        names=["name"],
    )
    return viewpoints["name"].dropna().astype(str).drop_duplicates().tolist()


def parquet_files(path: FilePath) -> list[pathlib.Path]:
    """Return parquet files represented by a file path or directory path."""
    path = pathlib.Path(path)
    if path.is_dir():
        return sorted(path.glob("*.parquet"))
    return [path]


def _normalise_nullable_viewpoints(reporters_df):
    import pandas as pd

    return reporters_df["viewpoint"].astype("string").replace(
        {"": pd.NA, "None": pd.NA, "nan": pd.NA, "<NA>": pd.NA}
    )


def _validate_reporter_columns(reporters_df, parquet_file: pathlib.Path) -> None:
    required_columns = {"viewpoint"}
    missing_columns = required_columns - set(reporters_df.columns)
    if missing_columns:
        raise ValueError(
            f"Reporter file {parquet_file} is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}"
        )


def write_countable_reporters(
    reporters: FilePath, viewpoint_path: FilePath, output_dir: FilePath
) -> pathlib.Path:
    """Write reporter parquet files with viewpoint categories limited to real baits.

    ``capcruncher-tools`` expects the reporter ``viewpoint`` category set to contain
    only viewpoints from the bait BED. Older CapCruncher reporter files can carry
    unused synthetic categories, so this normalises categories while still rejecting
    actual non-viewpoint values.
    """
    import pandas as pd

    valid_viewpoints = valid_viewpoint_names(viewpoint_path)
    if not valid_viewpoints:
        raise ValueError(f"No viewpoints found in {viewpoint_path}")

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for index, parquet_file in enumerate(parquet_files(reporters)):
        reporters_df = pd.read_parquet(parquet_file)
        _validate_reporter_columns(reporters_df, parquet_file)
        reporters_df["viewpoint"] = _normalise_nullable_viewpoints(reporters_df)
        invalid_viewpoints = sorted(
            set(reporters_df["viewpoint"].dropna().astype(str)) - set(valid_viewpoints)
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
    """Count reporter interactions using the external ``capcruncher-tools`` API."""
    from capcruncher_tools.api import count_interactions as count_interactions_records

    options = InteractionCountOptions(
        reporters=reporters,
        output=output,
        remove_exclusions=remove_exclusions,
        remove_viewpoint=remove_viewpoint,
        subsample=subsample,
        fragment_map=fragment_map,
        viewpoint_path=viewpoint_path,
        n_cores=n_cores,
        assay=assay,
        executor=executor,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        countable_reporters = write_countable_reporters(
            reporters=options.reporters,
            viewpoint_path=options.viewpoint_path,
            output_dir=tmpdir,
        )

        clr = count_interactions_records(
            reporters=countable_reporters,
            output=options.output,
            remove_exclusions=options.remove_exclusions,
            remove_viewpoint=options.remove_viewpoint,
            subsample=options.subsample,
            fragment_map=options.fragment_map,
            viewpoint_path=options.viewpoint_path,
            n_cores=options.n_cores,
            assay=options.assay,
            executor=options.executor,
            **kwargs,
        )

    return clr
