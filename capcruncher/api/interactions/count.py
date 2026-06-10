import os
import tempfile
from pathlib import Path
from typing import Any, cast
from uuid import uuid4

from loguru import logger
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, field_validator
from tqdm import tqdm

from capcruncher.api.interactions.pixels import iter_count_results
from capcruncher.api.interactions.reporters import (
    summarise_reporter_viewpoints,
    write_countable_reporters,
)
from capcruncher.types import (
    VALID_ASSAYS,
    VALID_EXECUTORS,
    Assay,
    Executor,
    validate_choice,
)


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

        reporters_for_counting = Path(countable_reporters)
        reporter_summary = summarise_reporter_viewpoints(reporters_for_counting)

        if options.fragment_map is None:
            raise ValueError("fragment_map is required.")

        import polars as pl

        bins = pl.read_csv(
            options.fragment_map,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "name"],
            schema_overrides={"chrom": pl.String},
        ).to_pandas()
        bins["chrom"] = bins["chrom"].astype("category")

        count_kwargs = {
            "parquet": os.fspath(reporters_for_counting / "*.parquet"),
            "remove_exclusions": options.remove_exclusions,
            "remove_viewpoint": options.remove_viewpoint,
            "subsample": options.subsample,
            "low_memory": reporter_summary.low_memory,
            "partitions": reporter_summary.partitions,
        }

        counted_viewpoints = []
        for viewpoint in reporter_summary.viewpoints:
            if reporter_summary.viewpoint_sizes.get(viewpoint, 0) == 0:
                logger.warning(
                    "No reporter rows found for viewpoint category "
                    f"{viewpoint}; skipping cooler creation"
                )
            else:
                counted_viewpoints.append(viewpoint)

        coolers = []
        for viewpoint, counts in tqdm(
            iter_count_results(
                counted_viewpoints,
                count_kwargs,
                executor_value,
                options.n_cores,
            ),
            total=len(counted_viewpoints),
        ):
            if counts.empty:
                logger.warning(
                    "No interactions found for viewpoint "
                    f"{viewpoint}; skipping cooler creation"
                )
                continue

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
        if not coolers:
            logger.warning(
                "No non-empty interaction counts were generated; "
                f"writing an empty HDF5 container at {options.output}"
            )
        merge_coolers(coolers, output=Path(options.output))

    return os.fspath(options.output)
