import cooler
import ujson as json
import pathlib


def viewpoints_with_interactions(cooler_path):
    viewpoints = []
    for viewpoint in cooler.api.list_coolers(cooler_path):
        clr = cooler.Cooler(f"{cooler_path}::{viewpoint}")
        count = clr.pixels()[:]["count"].sum()
        if count > 0:
            viewpoints.append(viewpoint)
    return viewpoints


def write_viewpoints_with_interactions(coolers, samples, outdir):
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for sample, cooler_path in zip(samples, coolers):
        viewpoints = viewpoints_with_interactions(cooler_path)
        with open(outdir / f"{sample}.json", "w") as f:
            json.dump(viewpoints, f)


if "snakemake" in globals():
    write_viewpoints_with_interactions(
        globals()["snakemake"].input[0],
        globals()["snakemake"].params.samples,
        globals()["snakemake"].params.outdir,
    )
