import pathlib
import shutil


def repartition_filtered_slices(
    flashed_slices, pe_slices, output_flashed, output_pe, log_path
):
    output_flashed = pathlib.Path(output_flashed)
    output_pe = pathlib.Path(output_pe)
    log_path = pathlib.Path(log_path)

    output_flashed.mkdir(parents=True, exist_ok=True)
    output_pe.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with open(log_path, "w") as log_handle:
        for source in flashed_slices:
            log_handle.write(f"Moving {source} to {output_flashed}\n")
            shutil.move(source, output_flashed)
        for source in pe_slices:
            log_handle.write(f"Moving {source} to {output_pe}\n")
            shutil.move(source, output_pe)


def main(snakemake):
    repartition_filtered_slices(
        snakemake.input.flashed,
        snakemake.input.pe,
        snakemake.output.slices_flashed,
        snakemake.output.slices_pe,
        snakemake.log[0],
    )


if "snakemake" in globals():
    main(globals()["snakemake"])
