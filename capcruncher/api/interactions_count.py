import tempfile
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, PositiveInt, field_validator
from capcruncher.types import Assay, Executor


class InteractionCountOptions(BaseModel):
    """Validated options for reporter interaction counting."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    reporters: Path
    output: Path = Path("CC_cooler.hdf5")
    remove_exclusions: bool = False
    remove_viewpoint: bool = False
    subsample: float = Field(default=0, ge=0, le=1)
    fragment_map: Path | None = None
    viewpoint_path: Path
    n_cores: PositiveInt = 1
    assay: Assay = Assay.CAPTURE
    executor: Executor = Executor.LOCAL

    @field_validator("reporters", "fragment_map", "viewpoint_path")
    @classmethod
    def existing_input_path(cls, value: Path | None) -> Path | None:
        if value is not None and not value.exists():
            raise ValueError(f"Input path does not exist: {value}")
        return value


def valid_viewpoint_names(viewpoint_path: Path | str) -> list[str]:
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


def parquet_files(path: Path | str) -> list[Path]:
    """Return parquet files represented by a file path or directory path."""
    path = Path(path)
    if path.is_dir():
        return sorted(path.glob("*.parquet"))
    return [path]


def _normalise_nullable_viewpoints(reporters_df):
    import pandas as pd

    return reporters_df["viewpoint"].astype("string").replace(
        {"": pd.NA, "None": pd.NA, "nan": pd.NA, "<NA>": pd.NA}
    )


def _validate_reporter_columns(reporters_df, parquet_file: Path) -> None:
    required_columns = {"viewpoint"}
    missing_columns = required_columns - set(reporters_df.columns)
    if missing_columns:
        raise ValueError(
            f"Reporter file {parquet_file} is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}"
        )


def write_countable_reporters(
    reporters: Path | str, viewpoint_path: Path | str, output_dir: Path | str
) -> Path:
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

    output_dir = Path(output_dir)
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
    reporters: Path | str,
    output: Path | str = Path("CC_cooler.hdf5"),
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    fragment_map: Path | str | None = None,
    viewpoint_path: Path | str | None = None,
    n_cores: int = 1,
    assay: Assay = Assay.CAPTURE,
    executor: Executor = Executor.LOCAL,
    **kwargs: Any,
) -> Path | str:
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
            reporters=str(countable_reporters),
            output=str(options.output),
            remove_exclusions=options.remove_exclusions,
            remove_viewpoint=options.remove_viewpoint,
            subsample=options.subsample,
            fragment_map=str(options.fragment_map) if options.fragment_map else None,
            viewpoint_path=str(options.viewpoint_path),
            n_cores=options.n_cores,
            assay=options.assay.value,
            executor=options.executor.value,
            **kwargs,
        )

    return clr
