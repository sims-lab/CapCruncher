import os

import pandas as pd
import pytest
import subprocess
from click.testing import CliRunner
from capcruncher.cli import cli
from capcruncher.utils import (
    bed_has_duplicate_names,
    bed_has_name,
    convert_bed_to_pr,
    convert_bed_to_dataframe,
    convert_interval_to_coords,
    format_coordinates,
    intersect_bins,
    is_valid_bed,
)


@pytest.fixture(scope="module")
def cli_runner():
    return CliRunner()


@pytest.fixture(scope="module")
def data_path_pipeline():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "data_for_pipeline_run")
    return data_dir


@pytest.fixture(scope="module")
def data_path_utils():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "utils")
    return data_dir


@pytest.fixture(scope="module")
def data_path_alignment_annotation():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "alignment_annotation")
    return data_dir


@pytest.fixture(scope="module")
def genome(data_path_pipeline):
    return os.path.join(data_path_pipeline, "chr14.fa.gz")


@pytest.fixture(scope="module")
def indicies(data_path_pipeline, genome):

    indicies = os.path.join(data_path_pipeline, "chr14_bowtie2_indicies")
    if not os.path.exists(indicies):
        os.mkdir(indicies)
        cmd = f"bowtie2-build {genome} {indicies}/bt2 --threads 8"
        subprocess.run(cmd.split())

    return os.path.join(indicies, "bt2")


@pytest.mark.parametrize(
    "viewpoints,flags",
    [
        (
            "viewpoints.fa",
            ["-r", "dpnii"],
        ),
    ],
)
def test_viewpoint_coordinates(
    cli_runner,
    viewpoints,
    genome,
    indicies,
    data_path_pipeline,
    data_path_utils,
    flags,
    tmpdir,
):

    viewpoints = os.path.join(data_path_utils, viewpoints)
    outfile = os.path.join(tmpdir, "viewpoint_coords.bed")

    result = cli_runner.invoke(
        cli,
        [
            "utilities",
            "viewpoint-coordinates",
            "-v",
            viewpoints,
            "-i",
            indicies,
            "-g",
            genome,
            "-o",
            outfile,
            *flags,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)


def test_bed_validation_and_formatting(data_path_alignment_annotation):
    capture = os.path.join(data_path_alignment_annotation, "test_capture.bed")
    slices = os.path.join(data_path_alignment_annotation, "test_slices_sorted.bed")
    blank = os.path.join(data_path_alignment_annotation, "blank.bed")

    assert is_valid_bed(capture)
    assert bed_has_name(capture)
    assert not bed_has_duplicate_names(capture)
    assert not is_valid_bed(blank)
    assert bed_has_duplicate_names(
        pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [0, 10],
                "end": [5, 15],
                "name": ["dup", "dup"],
            }
        )
    )

    coords = format_coordinates("chr1:10-20")
    assert convert_bed_to_dataframe(coords).loc[0, "name"] == "region_0"

    capture_pr = convert_bed_to_pr(capture)
    assert convert_bed_to_dataframe(capture_pr).loc[0, "name"] == "CAPTURE"

    named = format_coordinates(slices)
    assert convert_bed_to_dataframe(named).shape[0] == 4


def test_intersect_bins_and_interval_conversion():
    left = pd.DataFrame(
        {"chrom": ["chr1"], "start": [10], "end": [30], "name": ["left"]}
    )
    right = pd.DataFrame(
        {"chrom": ["chr1"], "start": [20], "end": [40], "name": ["right"]}
    )

    intersection = intersect_bins(left, right)
    assert list(intersection.columns) == [
        "chrom_1",
        "start_1",
        "end_1",
        "name_1",
        "chrom_2",
        "start_2",
        "end_2",
        "name_2",
        "overlap",
    ]
    assert intersection.loc[0, "overlap"] == 10

    name, coord = convert_interval_to_coords({"chrom": "chr1", "start": 10, "end": 30})
    assert name == "Unnammed"
    assert coord == "chr1:10-30"
