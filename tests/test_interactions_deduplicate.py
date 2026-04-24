import json

import polars as pl

from capcruncher.cli.interactions_deduplicate import deduplicate


def test_deduplicate_flashed_accepts_categorical_coordinates(tmp_path):
    slices = tmp_path / "slices.parquet"
    output = tmp_path / "deduplicated"
    statistics = tmp_path / "deduplication.json"

    pl.DataFrame(
        {
            "slice_id": [1, 2, 3, 4],
            "parent_id": [10, 10, 20, 20],
            "coordinates": ["chr1:1-10", "chr1:20-30", "chr1:1-10", "chr1:20-30"],
        },
        schema_overrides={"coordinates": pl.Categorical},
    ).write_parquet(slices)

    deduplicate(
        slices=slices,
        output=output,
        read_type="flashed",
        sample_name="sample-a",
        statistics=statistics,
    )

    assert list(output.rglob("*.parquet"))
    stats = json.loads(statistics.read_text())
    assert {
        key: stats[key]
        for key in [
            "sample",
            "read_type",
            "n_total_reads",
            "n_unique_reads",
            "n_total_slices",
            "n_unique_slices",
        ]
    } == {
        "sample": "sample-a",
        "read_type": "flashed",
        "n_total_reads": 2,
        "n_unique_reads": 1,
        "n_total_slices": 4,
        "n_unique_slices": 2,
    }
