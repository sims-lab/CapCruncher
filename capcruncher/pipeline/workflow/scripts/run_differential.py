import pathlib
import subprocess


EXPECTED_EMPTY_MARKERS = (
    "No differential interactions found",
    "No objects to concatenate",
    "Found array with 0 sample",
    "Found array with 0 feature",
)


def run_differential(
    counts,
    design_matrix,
    output_prefix,
    viewpoint,
    contrast,
    viewpoint_distance,
    output_dir,
    log_path,
):
    output_dir = pathlib.Path(output_dir)
    log_path = pathlib.Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "capcruncher",
        "interactions",
        "differential",
        *[str(count) for count in counts],
        "--design-matrix",
        str(design_matrix),
        "-o",
        str(output_prefix),
        "-v",
        str(viewpoint),
        "-c",
        str(contrast),
        "--viewpoint-distance",
        str(viewpoint_distance),
    ]

    with open(log_path, "w") as log_handle:
        completed = subprocess.run(
            cmd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if completed.returncode == 0:
            output_dir.mkdir(parents=True, exist_ok=True)
            return

        log_handle.flush()

    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    if any(marker in log_text for marker in EXPECTED_EMPTY_MARKERS):
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a", encoding="utf-8") as log_handle:
            log_handle.write(f"\nNo differential interactions found for {viewpoint}\n")
        return

    raise subprocess.CalledProcessError(completed.returncode, cmd)


def main(snakemake):
    run_differential(
        snakemake.input.counts,
        snakemake.input.design_matrix,
        snakemake.params.output_prefix,
        snakemake.params.viewpoint,
        snakemake.params.contrast,
        snakemake.params.viewpoint_distance,
        snakemake.output[0],
        snakemake.log[0],
    )


if "snakemake" in globals():
    main(globals()["snakemake"])
