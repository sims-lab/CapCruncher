import pathlib
import shutil


def combine_deduplicated_slices(source_dir, output_slices, log_path):
    source_dir = pathlib.Path(source_dir)
    output_slices = pathlib.Path(output_slices)
    log_path = pathlib.Path(log_path)

    output_slices.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with open(log_path, "w") as log_handle:
        for combined in ("flashed", "pe"):
            for source in (source_dir / combined).glob("*.parquet"):
                destination = output_slices / f"{combined}-{source.name}"
                log_handle.write(f"Moving {source} to {destination}\n")
                shutil.move(source, destination)


def main(snakemake):
    combine_deduplicated_slices(
        snakemake.params.source_dir,
        snakemake.output.slices,
        snakemake.log[0],
    )


if "snakemake" in globals():
    main(globals()["snakemake"])
