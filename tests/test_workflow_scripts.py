import importlib.util
import sys
import types
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


def test_make_ucsc_hub_builds_tracknado_metadata(tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")
    viewpoints = tmp_path / "viewpoints.bigBed"

    records = script.build_track_metadata(
        bigwigs=[
            tmp_path / "raw" / "SAMPLE-A_REP1_Slc25A37.bigWig",
            tmp_path / "norm" / "SAMPLE-A_REP1_Slc25A37.bigWig",
        ],
        bigwigs_summary=[
            tmp_path / "SAMPLE-A.mean-summary.Slc25A37.bigWig",
        ],
        bigwigs_comparison=[
            tmp_path / "SAMPLE-A-SAMPLE-B.mean-subtraction.Slc25A37.bigWig",
        ],
        viewpoints=viewpoints,
    )

    assert [
        {
            "category": record["category"],
            "normalisation": record["normalisation"],
            "sample": record["sample"],
            "aggregation": record["aggregation"],
            "ext": record["ext"],
        }
        for record in records
    ] == [
        {
            "category": "Replicates",
            "normalisation": "raw",
            "sample": "SAMPLE-A_REP1",
            "aggregation": "replicate",
            "ext": "bigWig",
        },
        {
            "category": "Replicates",
            "normalisation": "norm",
            "sample": "SAMPLE-A_REP1",
            "aggregation": "replicate",
            "ext": "bigWig",
        },
        {
            "category": "Aggregated",
            "normalisation": "norm",
            "sample": "SAMPLE-A",
            "aggregation": "mean",
            "ext": "bigWig",
        },
        {
            "category": "Subtraction",
            "normalisation": "norm",
            "sample": "SAMPLE-A-SAMPLE-B",
            "aggregation": "mean",
            "ext": "bigWig",
        },
        {
            "category": "Annotation",
            "normalisation": "viewpoints",
            "sample": "viewpoints",
            "aggregation": "viewpoints",
            "ext": "bigBed",
        },
    ]
    assert "overlay" not in records[-1]


def test_make_ucsc_hub_rejects_unsupported_track_names(tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")

    with pytest.raises(ValueError, match="Could not parse CapCruncher track path"):
        script.capcruncher_track_metadata(tmp_path / "sample.invalid-name.bigWig")


def test_make_ucsc_hub_uses_modern_tracknado_builder(monkeypatch, tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")
    calls = []

    class DummyHub:
        pass

    class DummyBuilder:
        def add_tracks(self, tracks):
            calls.append(("add_tracks", tracks))
            return self

        def with_metadata_extractor(self, extractor):
            calls.append(("with_metadata_extractor", extractor.__name__))
            return self

        def group_by(self, *columns, as_supertrack=False):
            calls.append(("group_by", columns, as_supertrack))
            return self

        def overlay_by(self, *columns):
            calls.append(("overlay_by", columns))
            return self

        def color_by(self, column):
            calls.append(("color_by", column))
            return self

        def with_custom_genome(
            self, name, twobit_file, organism, default_position
        ):
            calls.append(
                ("with_custom_genome", name, twobit_file, organism, default_position)
            )
            return self

        def build(self, **kwargs):
            calls.append(("build", kwargs))
            return DummyHub()

    monkeypatch.setitem(
        sys.modules,
        "tracknado",
        types.SimpleNamespace(HubBuilder=DummyBuilder),
    )

    result = script.build_hub(
        tracks=script.collect_track_paths(
            bigwigs=[tmp_path / "raw" / "SAMPLE-A_REP1_Slc25A37.bigWig"],
            bigwigs_summary=[],
            bigwigs_comparison=[],
            viewpoints=tmp_path / "viewpoints.bigBed",
        ),
        color_by="sample",
        genome="mm10",
        hub_name="capcruncher",
        hub_email="test@example.org",
        custom_genome=True,
        genome_twobit=tmp_path / "genome.2bit",
        genome_organism="Mouse",
        genome_default_position="chr1:1-100",
        report=tmp_path / "report.html",
        outdir=tmp_path / "hub",
    )

    assert isinstance(result, DummyHub)
    assert calls[1:] == [
        ("with_metadata_extractor", "capcruncher_track_metadata"),
        ("group_by", ("category", "normalisation"), True),
        ("group_by", ("sample", "viewpoint", "aggregation"), False),
        ("overlay_by", ("overlay",)),
        ("color_by", "sample"),
        (
            "with_custom_genome",
            "mm10",
            tmp_path / "genome.2bit",
            "Mouse",
            "chr1:1-100",
        ),
        (
            "build",
            {
                "name": "capcruncher",
                "genome": "mm10",
                "outdir": tmp_path / "hub",
                "hub_email": "test@example.org",
                "description_html": tmp_path / "report.html",
            },
        ),
    ]


def test_make_ucsc_hub_requires_twobit_for_custom_genome(tmp_path):
    script = load_workflow_script("make_ucsc_hub.py")

    with pytest.raises(ValueError, match="genome twoBit file"):
        script.build_hub(
            tracks=[
                tmp_path / "raw" / "SAMPLE-A_REP1_Slc25A37.bigWig",
                tmp_path / "viewpoints.bigBed",
            ],
            color_by="sample",
            genome="custom",
            hub_name="capcruncher",
            hub_email="test@example.org",
            custom_genome=True,
            genome_twobit=None,
            report=tmp_path / "report.html",
            outdir=tmp_path / "hub",
        )
