import os
import pandas as pd
import subprocess
import pyarrow.dataset as ds
import pathlib
from loguru import logger


def remove_duplicate_coordinates(
    slices_directory,
    output_slices,
    output_statistics,
    read_type,
    sample_name,
    log_path,
):
    dataset = ds.dataset(slices_directory, format="parquet")
    n_rows = dataset.count_rows()

    if n_rows != 0:
        cmd = [
            "capcruncher",
            "interactions",
            "deduplicate",
            slices_directory,
            "-o",
            output_slices,
            "--read-type",
            read_type,
            "--sample-name",
            sample_name,
            "--statistics",
            output_statistics,
        ]

        with open(log_path, "w") as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=f)

    else:
        logger.warning("The input dataset is empty, skipping deduplication step.")

        outdir = pathlib.Path(output_slices)

        logger.warning(f"Creating empty output directory: {outdir}")
        outdir.mkdir(parents=True, exist_ok=True)

        logger.warning(f"Creating empty stats file: {output_statistics}")
        pd.DataFrame().to_csv(output_statistics)


def main(snakemake):
    try:
        remove_duplicate_coordinates(
            slices_directory=snakemake.input.slices_directory,
            output_slices=snakemake.output.slices,
            output_statistics=snakemake.output.statistics,
            read_type=snakemake.params.read_type,
            sample_name=snakemake.params.sample_name,
            log_path=snakemake.log[0],
        )
    except Exception as exc:
        print(exc)
        os.makedirs(snakemake.output.slices, exist_ok=True)


if "snakemake" in globals():
    main(globals()["snakemake"])
