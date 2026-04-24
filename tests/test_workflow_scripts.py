import importlib.util
from pathlib import Path

import polars as pl
import pytest


def load_workflow_script(script_name):
    script_path = (
        Path(__file__).resolve().parents[1]
        / "capcruncher"
        / "pipeline"
        / "workflow"
        / "scripts"
        / script_name
    )
    spec = importlib.util.spec_from_file_location(script_name, script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_validation_confirm_annotated_viewpoints_present_counts_current_polars(
    tmp_path,
):
    script = load_workflow_script("validation_confirm_annotated_viewpoints_present.py")
    slices_a = tmp_path / "slices_a.parquet"
    slices_b = tmp_path / "slices_b.parquet"
    viewpoints = tmp_path / "viewpoints.bed"
    counts = tmp_path / "viewpoints_present.tsv"
    sentinel = tmp_path / "validated.sentinel"

    pl.DataFrame({"capture": ["vp1", "vp1", "vp2"]}).write_parquet(slices_a)
    pl.DataFrame({"capture": ["vp2"]}).write_parquet(slices_b)
    viewpoints.write_text("chr1\t0\t10\tvp1\nchr1\t20\t30\tvp2\n")

    script.validate_viewpoints_present(
        [slices_a, slices_b],
        viewpoints,
        counts,
        sentinel,
    )

    df_counts = pl.read_csv(counts, separator="\t")
    assert sentinel.exists()
    assert df_counts.select(["capture", "n_slices"]).sort("capture").to_dicts() == [
        {"capture": "vp1", "n_slices": 2},
        {"capture": "vp2", "n_slices": 2},
    ]


def test_validation_confirm_annotated_viewpoints_present_reports_missing(tmp_path):
    script = load_workflow_script("validation_confirm_annotated_viewpoints_present.py")
    slices = tmp_path / "slices.parquet"
    viewpoints = tmp_path / "viewpoints.bed"

    pl.DataFrame({"capture": ["vp1"]}).write_parquet(slices)
    viewpoints.write_text("chr1\t0\t10\tvp1\nchr1\t20\t30\tvp2\n")

    with pytest.raises(ValueError, match="vp2"):
        script.validate_viewpoints_present(
            [slices],
            viewpoints,
            tmp_path / "viewpoints_present.tsv",
            tmp_path / "validated.sentinel",
        )


def test_count_identified_viewpoints_filters_empty_values(tmp_path):
    script = load_workflow_script("count_identified_viewpoints.py")
    slices_dir = tmp_path / "slices"
    slices_dir.mkdir()
    output = tmp_path / "identified_viewpoints.tsv"

    pl.DataFrame(
        {
            "viewpoint": ["vp1", "vp1", "vp2", "", None],
            "pe": ["pe1", "pe1", "", "pe2", "pe3"],
        }
    ).write_parquet(slices_dir / "slices.parquet")

    script.write_identified_viewpoints(slices_dir, output)

    assert pl.read_csv(output, separator="\t").sort("viewpoint").to_dicts() == [
        {"viewpoint": "vp1", "pe": "pe1"},
    ]
