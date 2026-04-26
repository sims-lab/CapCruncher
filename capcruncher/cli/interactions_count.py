import os
import pathlib
import tempfile
from typing import Literal

import pandas as pd

type FilePath = str | os.PathLike[str]


def _valid_viewpoint_names(viewpoint_path: FilePath) -> list[str]:
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


def count(
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
    **kwargs,
) -> FilePath:
    """
    Counts interactions between the viewpoint and the rest of the genome.

    Args:
        reporters: Path to reporters file.
        output: Output file name.
        remove_exclusions: Remove excluded regions.
        remove_viewpoint: Remove capture regions.
        subsample: Subsample reads.
        fragment_map: Path to fragment map.
        viewpoint_path: Path to viewpoint file.
        n_cores: Number of cores.
        assay: Assay type.
        executor: Runtime used for per-viewpoint counting.
        **kwargs: Additional arguments.
    Returns:
        Path to the generated cooler file.

    """
    from capcruncher_tools.api import count_interactions

    with tempfile.TemporaryDirectory() as tmpdir:
        countable_reporters = _write_countable_reporters(
            reporters=reporters,
            viewpoint_path=viewpoint_path,
            output_dir=tmpdir,
        )

        clr = count_interactions(
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
