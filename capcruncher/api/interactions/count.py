import os
import tempfile
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path
from typing import Any, cast
from uuid import uuid4

import pandas as pd
from loguru import logger
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, field_validator
from tqdm import tqdm

from capcruncher.types import Assay, Executor, VALID_ASSAYS, VALID_EXECUTORS, validate_choice


class InteractionCountOptions(BaseModel):
    """Validated options for reporter interaction counting."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    reporters: Path | str
    output: Path | str = Path("CC_cooler.hdf5")
    remove_exclusions: bool = False
    remove_viewpoint: bool = False
    subsample: float = Field(default=0, ge=0, le=1)
    fragment_map: Path | str | None = None
    viewpoint_path: Path | str
    n_cores: PositiveInt = 1
    assay: Assay | str = Assay.CAPTURE
    executor: Executor | str = Executor.LOCAL

    @field_validator("reporters", "fragment_map", "viewpoint_path", mode="before")
    @classmethod
    def existing_input_path(cls, value: Path | str | None) -> Path | None:
        if value is None:
            return None
        path = Path(value)
        if not path.exists():
            raise ValueError(f"Input path does not exist: {path}")
        return path

    @field_validator("assay", mode="before")
    @classmethod
    def validate_assay(cls, value: Assay | str) -> Assay:
        return validate_choice(value, VALID_ASSAYS, "assay")

    @field_validator("executor", mode="before")
    @classmethod
    def validate_executor(cls, value: Executor | str) -> Executor:
        return validate_choice(value, VALID_EXECUTORS, "executor")


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
    assay: Assay | str = Assay.CAPTURE,
    executor: Executor | str = Executor.LOCAL,
    **kwargs: Any,
) -> Path | str:
    """Count reporter interactions and write CapCruncher cooler output."""
    import pyranges1 as pr

    from capcruncher.api.interactions.cooler.create import create_cooler_cc
    from capcruncher.api.interactions.cooler.merge import merge_coolers

    if viewpoint_path is None:
        raise ValueError("viewpoint_path is required.")

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

        assay_value = cast(Assay, options.assay).value
        executor_value = cast(Executor, options.executor).value

        logger.info("Extracting viewpoint names and sizes")
        reporters_for_counting = Path(countable_reporters)
        viewpoint_df = pd.read_parquet(
            reporters_for_counting, engine="pyarrow", columns=["viewpoint"]
        )
        if hasattr(viewpoint_df["viewpoint"], "cat") and hasattr(
            viewpoint_df["viewpoint"].cat, "categories"
        ):
            viewpoints = viewpoint_df["viewpoint"].cat.categories.to_list()
        else:
            viewpoints = (
                viewpoint_df["viewpoint"].dropna().drop_duplicates().astype(str).to_list()
            )
        viewpoint_sizes = viewpoint_df["viewpoint"].value_counts()
        viewpoint_sizes_dict = viewpoint_sizes.to_dict()
        viewpoint_sizes_df = pd.DataFrame.from_dict(
            viewpoint_sizes_dict, orient="index", columns=["n_slices"]
        )

        logger.info(f"Number of viewpoints: {len(viewpoints)}")
        logger.info(f"Number of slices per viewpoint:\n{viewpoint_sizes_df}")

        if any(size > 2e6 for size in viewpoint_sizes_dict.values()):
            logger.warning(
                "High number of slices per viewpoint detected. Switching to low memory mode"
            )
            low_memory = True
            partition_df = pd.read_parquet(
                reporters_for_counting, engine="pyarrow", columns=["bam"]
            )
            partitions = partition_df["bam"].cat.categories.to_list()
        else:
            low_memory = False
            partitions = None

        if options.fragment_map is None:
            raise ValueError("fragment_map is required.")

        bins = pr.read_bed(Path(options.fragment_map)).rename(
            columns={
                "Chromosome": "chrom",
                "Start": "start",
                "End": "end",
                "Name": "name",
            }
        )
        bins["chrom"] = bins["chrom"].astype("string").astype("category")

        count_kwargs = {
            "parquet": os.fspath(reporters_for_counting / "*.parquet"),
            "remove_exclusions": options.remove_exclusions,
            "remove_viewpoint": options.remove_viewpoint,
            "subsample": options.subsample,
            "low_memory": low_memory,
            "partitions": partitions,
        }

        coolers = []
        for viewpoint, counts in tqdm(
            _iter_count_results(
                viewpoints, count_kwargs, executor_value, options.n_cores
            ),
            total=len(viewpoints),
        ):
            cooler_uri = create_cooler_cc(
                output_prefix=os.fspath(Path(tmpdir) / f"{uuid4().hex}.hdf5"),
                pixels=counts,
                bins=bins,
                viewpoint_name=viewpoint,
                viewpoint_path=Path(options.viewpoint_path),
                assay=assay_value,
                **kwargs,
            )
            coolers.append(cooler_uri.split("::")[0])

        logger.info(f"Making final cooler at {options.output}")
        merge_coolers(coolers, output=Path(options.output))

    return os.fspath(options.output)


def _iter_count_results(
    viewpoints: Iterable[str],
    count_kwargs: dict[str, Any],
    executor: str,
    n_cores: int,
) -> Iterable[tuple[str, pd.DataFrame]]:
    from capcruncher_tools.count import count_viewpoint_pixels

    if executor == Executor.LOCAL.value:
        for viewpoint in viewpoints:
            yield count_viewpoint_pixels(viewpoint=viewpoint, **count_kwargs)
        return

    if executor == Executor.PROCESS.value:
        process_kwargs: dict[str, Any] = {"max_workers": n_cores}
        try:
            process_kwargs["mp_context"] = get_context("fork")
        except ValueError:
            pass
        try:
            with ProcessPoolExecutor(**process_kwargs) as pool:
                futures = [
                    pool.submit(
                        count_viewpoint_pixels, viewpoint=viewpoint, **count_kwargs
                    )
                    for viewpoint in viewpoints
                ]
                for future in as_completed(futures):
                    yield future.result()
        except PermissionError:
            logger.warning(
                "Process executor is unavailable in this environment; falling back to local execution"
            )
            for viewpoint in viewpoints:
                yield count_viewpoint_pixels(viewpoint=viewpoint, **count_kwargs)
        return

    if executor == Executor.RAY.value:
        try:
            import ray
        except ImportError as exc:
            raise RuntimeError(
                "Ray executor requested but ray is not installed. Install capcruncher-tools[ray]."
            ) from exc

        ray.init(num_cpus=n_cores, ignore_reinit_error=True)
        remote_count = ray.remote(count_viewpoint_pixels)
        futures = [
            remote_count.remote(viewpoint=viewpoint, **count_kwargs)
            for viewpoint in viewpoints
        ]
        while futures:
            completed, futures = ray.wait(futures)
            for future in completed:
                yield ray.get(future)
        return

    raise ValueError(f"Unknown executor: {executor}")
