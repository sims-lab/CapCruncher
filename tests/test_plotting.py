import os
import pathlib
from typing import Any, Dict, List
from unittest import TestCase

import pytest

from capcruncher.api.plotting import (
    CCFigure,
    CCTrack,
)


def can_import_coolbox():
    try:
        import coolbox.api as cb

        return True
    except ImportError:
        return False


@pytest.fixture(scope="module")
def repo_path():
    fn = pathlib.Path(__file__).resolve()
    dirname = fn.parent
    return dirname.parent


@pytest.fixture(scope="module")
def package_path(repo_path):
    return repo_path.joinpath("capcruncher")


@pytest.fixture(scope="module")
def data_path():
    fn = pathlib.Path(__file__).resolve()
    dirname = fn.parent
    data_dir = dirname.joinpath("data")
    return data_dir


# Fixture to create the CCMatrix object for testing
@pytest.fixture(params=["raw", "n_interactions", "icen_cis", "icen_scale"])
def heatmap(data_path):
    file_path = data_path / "reporters_store" / "SAMPLE-A_REP1_binned.hdf5"
    track = CCTrack(
        file=str(file_path), binsize=5000, viewpoint="Slc25A37", file_type="heatmap"
    )
    return track


@pytest.fixture()
def heatmap_summary(data_path):
    file_path_1 = data_path / "reporters_store" / "SAMPLE-A_REP1_binned.hdf5"
    file_path_2 = data_path / "reporters_store" / "SAMPLE-A_REP1_binned.hdf5"
    track = CCTrack(
        file=[file_path_1, file_path_2],
        binsize=5000,
        viewpoint="Slc25A37",
        file_type="heatmap_summary",
    )
    return track


# Fixture to create the CCBigWig object for testing
@pytest.fixture
def bigwig(data_path):
    file_path = (
        data_path / "test_bigwigs" / "Slc25A37-test-1x_1.normalised.Slc25A37.bigWig"
    )
    track = CCTrack(file=str(file_path), file_type="bigwig")
    return track


@pytest.fixture
def bed(data_path):
    file_path = (
        data_path / "data_for_pipeline_run" / "mm9_capture_viewpoints_Slc25A37.bed"
    )
    return CCTrack(file=str(file_path), file_type="bed")


@pytest.fixture
def bigwig_summary(data_path):
    file_paths = (data_path / "test_bigwigs").glob("*1x*.bigWig")
    track = CCTrack(file=file_paths, file_type="bigwig_summary")
    return track


@pytest.fixture
def arcs(data_path):
    file_path = data_path / "plotting" / "test.bedpe"
    return CCTrack(file=str(file_path), file_type="Arcs")


@pytest.fixture
def coordinates():
    chrom = "chr14"
    start = 69878303
    end = 69946880
    return f"{chrom}:{start - 1e4: .0f}-{end + 1e4: .0f}"


@pytest.mark.skipif(can_import_coolbox() is False, reason="Coolbox not installed")
def test_plotting(
    tmpdir, heatmap, heatmap_summary, bigwig, bigwig_summary, bed, coordinates, arcs
):
    # Create the figure
    fig = CCFigure()
    # Add the matrix
    fig.add_track(heatmap)
    # Add the matrix collection
    fig.add_track(heatmap_summary)
    # Add the bigwig
    fig.add_track(bigwig)
    # Add the bigwig collection
    fig.add_track(bigwig_summary)
    # Add the scale bar
    fig.add_track(CCTrack(None, file_type="scale"))
    # Add the bed file
    fig.add_track(bed)
    # Add the x-axis
    fig.add_track(CCTrack(None, file_type="xaxis"))
    # Add a random coolbox track
    fig.add_track(arcs)

    # Save the figure
    fig.save(coordinates, output=tmpdir / "test_plotting.png")

    # Check the file exists
    assert (tmpdir / "test_plotting.png").exists()


def test_toml_conversion(tmpdir, bigwig):
    fig = CCFigure()
    bigwig.properties["title"] = "test"
    fig.add_track(bigwig)

    toml_path = tmpdir / "test_toml_conversion.toml"
    fig.to_toml(toml_path)

    assert toml_path.exists()

    fig2 = CCFigure.from_toml(toml_path)

    track = next(iter(fig2.tracks))
    assert track.file == bigwig.file
