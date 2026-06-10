import pathlib


def save_design(design, output, log_path):
    design.to_csv(output, sep="\t", index=False)
    log_path = pathlib.Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(f"Wrote design matrix to {output}\n", encoding="utf-8")


def main(snakemake):
    save_design(snakemake.params.design, snakemake.output[0], snakemake.log[0])


if "snakemake" in globals():
    main(globals()["snakemake"])
