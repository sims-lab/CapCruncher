import pytest
import os
import coolbox.api as cb
from coolbox.core.track import Track
from typing import List, Dict, Any
import pathlib
from capcruncher.api.plotting import (
    CCBigWig,
    CCMatrix,
    CCBigWigCollection,
    ScaleBar,
    CCSimpleBed,
    CCXAxisGenomic,
    CCTrack,
    CCFigure,
)


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
@pytest.fixture
def heatmap(data_path):
    file_path = data_path / "reporters_store" / "SAMPLE-A_REP1_binned.hdf5"
    track = CCTrack(
        file=str(file_path), binsize=5000, viewpoint="Slc25A37", file_type="heatmap"
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
def coordinates():
    return "chr14:69902454-69903469"


def test_plotting(tmpdir, heatmap, bigwig, bigwig_summary, bed, coordinates):
    # Create the figure
    fig = CCFigure()
    # Add the matrix
    fig.add_track(heatmap)
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

    # Save the figure
    fig.save(coordinates, output=tmpdir / "test_plotting.png")

    # Check the file exists
    assert (tmpdir / "test_plotting.png").exists()

    # Check the file is not empty
    assert (tmpdir / "test_plotting.png").stat().st_size > 0


def test_toml_conversion(tmpdir, bigwig):
    fig = CCFigure()
    fig.add_track(bigwig)

    toml_path = tmpdir / "test_toml_conversion.toml"
    fig.to_toml(toml_path)

    assert toml_path.exists()

    fig2 = CCFigure.from_toml(toml_path)

    assert fig2.tracks[0].file == bigwig.file
    assert fig2.tracks[0].name == bigwig.name
    assert fig2.tracks[0].type == bigwig.type
