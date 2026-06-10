"""Count re-ligation events and cis interaction distances from filtered slices.

Produces a JSON file with two keys:
  religation    – per-sample/viewpoint/read_type re-ligation counts
  cis_distances – per-sample/viewpoint/read_type cis distance histogram
"""

import json
import pathlib

import polars as pl
from loguru import logger

DISTANCE_BINS = [0, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, float("inf")]
BIN_LABELS = ["<1kb", "1kb-10kb", "10kb-100kb", "100kb-1Mb", "1Mb-10Mb", ">10Mb"]


def _read_type_label(value: str) -> str:
    return {"flashed": "Combined", "pe": "Non-Combined"}.get(value, value)


def _assign_distance_bin(distance: int) -> str:
    for i, upper in enumerate(DISTANCE_BINS[1:]):
        if distance < upper:
            return BIN_LABELS[i]
    return BIN_LABELS[-1]


def compute_stats(
    slices: pl.DataFrame,
    sample_name: str,
) -> tuple[list[dict], list[dict]]:
    # Viewpoint slices carry the capture site coordinates
    captures = (
        slices.filter(pl.col("capture_count") > 0)
        .group_by("parent_read", "viewpoint", "pe")
        .agg(
            pl.col("restriction_fragment").first().alias("capture_fragment"),
            pl.col("chrom").first().alias("capture_chrom"),
            pl.col("start").first().alias("capture_start"),
            pl.col("end").first().alias("capture_end"),
        )
    )

    reporters = slices.filter(pl.col("capture_count") == 0)

    if reporters.is_empty() or captures.is_empty():
        return [], []

    joined = reporters.join(captures, on=["parent_read", "pe"], how="left")

    religation_rows = []
    distance_rows = []

    viewpoints = captures["viewpoint"].drop_nulls().unique().to_list()

    for vp in sorted(viewpoints):
        vp_data = joined.filter(pl.col("viewpoint") == vp)
        if vp_data.is_empty():
            continue

        for read_type_raw in vp_data["pe"].unique().to_list():
            rt_data = vp_data.filter(pl.col("pe") == read_type_raw)
            read_type = _read_type_label(read_type_raw)

            n_total = rt_data.height
            n_relig = rt_data.filter(
                (pl.col("restriction_fragment") - pl.col("capture_fragment")).abs() == 1
            ).height

            religation_rows.append({
                "sample": sample_name,
                "viewpoint": vp,
                "read_type": read_type,
                "n_total_reporters": n_total,
                "n_religation": n_relig,
                "percentage_religation": round(n_relig / n_total * 100, 2) if n_total > 0 else 0.0,
            })

            # Cis distance distribution
            cis = rt_data.filter(pl.col("chrom") == pl.col("capture_chrom"))
            if cis.is_empty():
                continue

            cis_with_dist = cis.with_columns(
                (
                    ((pl.col("start") + pl.col("end")) / 2
                     - (pl.col("capture_start") + pl.col("capture_end")) / 2).abs()
                ).cast(pl.Int64).alias("distance")
            )

            bin_counts: dict[str, int] = {label: 0 for label in BIN_LABELS}
            for row in cis_with_dist["distance"].to_list():
                bin_counts[_assign_distance_bin(row)] += 1

            for i, (label, cnt) in enumerate(zip(BIN_LABELS, [bin_counts[b] for b in BIN_LABELS])):
                distance_rows.append({
                    "sample": sample_name,
                    "viewpoint": vp,
                    "read_type": read_type,
                    "distance_bin": label,
                    "distance_count": cnt,
                    "bin_order": i,
                })

    return religation_rows, distance_rows


def main(snakemake) -> None:
    sample_name: str = snakemake.params.sample_name
    slices_path: str = snakemake.input.slices
    output_path: str = snakemake.output.stats

    logger.info(f"Loading slices for {sample_name} from {slices_path}")
    slices = pl.read_parquet(slices_path)
    logger.info(f"Loaded {slices.height} slices")

    religation_rows, distance_rows = compute_stats(slices, sample_name)

    logger.info(
        f"Found {sum(r['n_religation'] for r in religation_rows)} re-ligation events "
        f"across {len(religation_rows)} viewpoint/read-type combinations"
    )

    out = pathlib.Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps({"religation": religation_rows, "cis_distances": distance_rows}),
        encoding="utf-8",
    )
    logger.info(f"Written stats to {output_path}")


if "snakemake" in globals():
    main(globals()["snakemake"])
