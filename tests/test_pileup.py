import pathlib

import pandas as pd
import pyranges1 as pr
import pytest

from capcruncher.api.interactions.bedgraph import CCBedgraph, CoolerBedGraph


@pytest.fixture(scope="module")
def data_dir():
    file = pathlib.Path(__file__).resolve()
    dirname = file.parent
    return dirname.joinpath("data", "reporters_store")


@pytest.fixture(scope="module")
def viewpoint():
    return "Slc25A37"


@pytest.fixture(scope="module")
def norm_regions():
    file = pathlib.Path(__file__).resolve()
    dirname = file.parent
    regions = dirname.joinpath("data", "data_for_pipeline_run", "regions_for_norm.bed")
    return regions


@pytest.fixture(scope="module")
def cooler_bedgraph(data_dir, viewpoint):
    # Create a CoolerBedGraph instance for testing
    uri = str(data_dir / f"SAMPLE-A_REP1.hdf5::/{viewpoint}")
    return CoolerBedGraph(uri)


def test_get_reporters(cooler_bedgraph):
    # Test the _get_reporters method
    reporters = cooler_bedgraph._get_reporters()
    assert isinstance(reporters, pd.DataFrame)
    assert "capture" in reporters.columns
    assert "reporter" in reporters.columns
    assert "count" in reporters.columns


@pytest.mark.parametrize("norm", ["raw", "n_cis"])
def test_extract_bedgraph(cooler_bedgraph, norm):
    # Test the extract_bedgraph method
    bedgraph = cooler_bedgraph.extract_bedgraph(normalisation=norm)
    assert isinstance(bedgraph, pd.DataFrame)
    assert "chrom" in bedgraph.columns
    assert "start" in bedgraph.columns
    assert "end" in bedgraph.columns
    assert "count" in bedgraph.columns


def test_norm_by_region(cooler_bedgraph, norm_regions):
    # Test the norm_by_region method
    bedgraph = cooler_bedgraph.extract_bedgraph(
        normalisation="region", region=norm_regions
    )
    assert isinstance(bedgraph, pd.DataFrame)
    assert "chrom" in bedgraph.columns
    assert "start" in bedgraph.columns
    assert "end" in bedgraph.columns


def test_to_pyranges(cooler_bedgraph):
    # Test the to_pyranges method
    pyranges = cooler_bedgraph.to_pyranges(normalisation="raw")
    assert isinstance(pyranges, pr.PyRanges)


def test_to_pyranges_expands_normalisation_kwargs(monkeypatch):
    bedgraph = CoolerBedGraph.__new__(CoolerBedGraph)
    captured = {}

    def fake_extract_bedgraph(self, normalisation="raw", **norm_kwargs):
        captured["normalisation"] = normalisation
        captured["norm_kwargs"] = norm_kwargs
        return pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [10],
                "end": [20],
                "count": [1.0],
            }
        )

    monkeypatch.setattr(CoolerBedGraph, "extract_bedgraph", fake_extract_bedgraph)

    converted = bedgraph.to_pyranges(normalisation="region", region="regions.bed")

    assert isinstance(converted, pr.PyRanges)
    assert captured == {
        "normalisation": "region",
        "norm_kwargs": {"region": "regions.bed"},
    }


def test_region_normalisation_uses_interval_overlap(tmp_path):
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t10\t20\tSlc25A37_region\n")

    bedgraph = CoolerBedGraph.__new__(CoolerBedGraph)
    bedgraph.viewpoint_name = "Slc25A37"
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [5, 30],
            "end": [15, 40],
            "count": [10.0, 5.0],
        }
    )

    bedgraph._normalise_by_regions(df, scale_factor=1e6, regions=regions)

    assert df["count"].tolist() == [1e6, 5e5]


def test_ccbedgraph_to_pyranges():
    bedgraph = CCBedgraph(
        df=pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [10],
                "end": [20],
                "score": [1.0],
            }
        )
    )

    converted = bedgraph.to_pyranges()

    assert isinstance(converted, pr.PyRanges)
