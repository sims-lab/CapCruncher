import pandas as pd
import pytest

import capcruncher.api.interactions.differential as differential
from capcruncher.api.interactions.differential import get_differential_interactions

try:
    import pydeseq2  # noqa: F401

    HAS_PYDESEQ2 = True
except ImportError:
    HAS_PYDESEQ2 = False


def test_get_differential_interactions_uses_pydeseq2_results_df(monkeypatch):
    captured = {}

    class FakeInference:
        def __init__(self, n_cpus):
            captured["n_cpus"] = n_cpus

    class FakeDeseqDataSet:
        def __init__(self, counts, metadata, design, refit_cooks, inference):
            captured["counts"] = counts
            captured["metadata"] = metadata
            captured["design"] = design
            captured["refit_cooks"] = refit_cooks
            captured["inference"] = inference
            self.obsm = {
                "design_matrix": pd.DataFrame(columns=["Intercept", "condition[T.B]"])
            }

        def deseq2(self):
            captured["deseq2"] = True

    class FakeDeseqStats:
        def __init__(self, dds, contrast, inference):
            captured["contrast"] = contrast
            captured["stats_inference"] = inference
            self.results_df = None

        def summary(self):
            captured["summary"] = True
            self.results_df = pd.DataFrame(
                {
                    "baseMean": [10.0],
                    "log2FoldChange": [1.25],
                    "lfcSE": [0.1],
                    "stat": [2.0],
                    "pvalue": [0.01],
                    "padj": [0.02],
                },
                index=["chr1:10-20"],
            )

        def lfc_shrink(self, coeff):
            captured["lfc_shrink_coeff"] = coeff
            self.results_df["log2FoldChange"] = 0.75

    monkeypatch.setattr(
        differential,
        "_load_pydeseq2",
        lambda: (FakeDeseqDataSet, FakeInference, FakeDeseqStats),
    )

    counts = pd.DataFrame(
        {"SAMPLE-A": [1.2], "SAMPLE-B": [4.8]},
        index=["chr1:10-20"],
    )
    design = pd.DataFrame(
        {"condition": ["A", "B"]},
        index=["SAMPLE-A", "SAMPLE-B"],
    )

    results = differential.get_differential_interactions(
        counts,
        design,
        contrast="condition",
        group_a="A",
        group_b="B",
        lfc_shrink=True,
    )

    assert captured["counts"].to_dict() == {
        "chr1:10-20": {"SAMPLE-A": 1.2, "SAMPLE-B": 4.8}
    }
    assert captured["design"] == "~condition"
    assert captured["contrast"] == ["condition", "B", "A"]
    assert captured["lfc_shrink_coeff"] == "condition[T.B]"
    assert results.loc["chr1:10-20", "log2FoldChange"] == 0.75
    assert results.loc["chr1:10-20", ["chrom", "start", "end"]].to_dict() == {
        "chrom": "chr1",
        "start": 10,
        "end": 20,
    }


def test_differential_reports_empty_counts_after_threshold(monkeypatch, tmp_path):
    monkeypatch.setattr(
        differential,
        "cooler_to_bedgraph",
        lambda **kwargs: pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [10],
                "end": [20],
                "count": [1],
            }
        ),
    )

    design = tmp_path / "design.tsv"
    pd.DataFrame({"sample": ["sample-a", "sample-b"], "condition": ["A", "B"]}).to_csv(
        design, sep="\t", index=False
    )

    with pytest.raises(ValueError, match="No differential interactions found"):
        differential.differential(
            interaction_files=["sample-a.hdf5", "sample-b.hdf5"],
            viewpoint="vp1",
            design_matrix=design,
            output_prefix=tmp_path / "differential" / "vp1",
            contrast="condition",
            viewpoint_distance=1000,
            threshold_count=20,
        )


@pytest.mark.skipif(not HAS_PYDESEQ2, reason="pydeseq2 not installed")
def test_differential_interactions_end_to_end():
    counts = pd.DataFrame(
        {
            "SAMPLE-A_REP1": [100, 5],
            "SAMPLE-A_REP2": [120, 8],
            "SAMPLE-B_REP1": [10, 200],
            "SAMPLE-B_REP2": [8, 180],
        },
        index=["chr1:100-200", "chr1:300-400"],
    )
    design = pd.DataFrame(
        {"condition": ["A", "A", "B", "B"]},
        index=["SAMPLE-A_REP1", "SAMPLE-A_REP2", "SAMPLE-B_REP1", "SAMPLE-B_REP2"],
    )

    results = get_differential_interactions(
        counts,
        design,
        contrast="condition",
        group_a="A",
        group_b="B",
        threshold_q=1.0,
    )

    assert isinstance(results, pd.DataFrame)
    for col in (
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "chrom",
        "start",
        "end",
    ):
        assert col in results.columns, f"Missing column: {col}"
