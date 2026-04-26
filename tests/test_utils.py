import os
import builtins
import importlib
import sys

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
    read_dataframes,
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


def test_interval_helpers_import_without_ray(monkeypatch):
    real_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == "ray" or name.startswith("ray."):
            raise ModuleNotFoundError("No module named 'ray'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    for module in (
        "capcruncher.api.pileup",
        "capcruncher.api.storage",
        "capcruncher.utils",
    ):
        importlib.reload(importlib.import_module(module))


def test_api_package_import_does_not_import_submodules(monkeypatch):
    for module in list(sys.modules):
        if module == "capcruncher.api" or module.startswith("capcruncher.api."):
            monkeypatch.delitem(sys.modules, module, raising=False)

    api = importlib.import_module("capcruncher.api")

    assert api.__name__ == "capcruncher.api"
    assert not any(
        module.startswith("capcruncher.api.") for module in sys.modules
    )


def test_differential_module_imports_without_pydeseq2(monkeypatch):
    real_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == "pydeseq2" or name.startswith("pydeseq2."):
            raise ModuleNotFoundError("No module named 'pydeseq2'", name="pydeseq2")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    for module in list(sys.modules):
        if module == "capcruncher.api.interactions_differential" or module.startswith(
            "pydeseq2"
        ):
            monkeypatch.delitem(sys.modules, module, raising=False)

    differential = importlib.import_module("capcruncher.api.interactions_differential")

    with pytest.raises(ModuleNotFoundError, match="differential"):
        differential._load_pydeseq2()


def test_read_dataframes_skips_empty_files(tmp_path):
    empty = tmp_path / "empty.tsv"
    nonempty = tmp_path / "nonempty.tsv"
    empty.touch()
    nonempty.write_text("sample\tvalue\nA\t1\n", encoding="utf-8")

    frames = read_dataframes([empty, nonempty], sep="\t")

    assert len(frames) == 1
    assert frames[0].to_dict("records") == [{"sample": "A", "value": 1}]


def test_read_dataframes_reports_all_empty_files(tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.touch()

    with pytest.raises(RuntimeError, match="All dataframes supplied are empty"):
        read_dataframes([empty], sep="\t")


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


def test_intersect_bins_no_overlap_returns_empty_schema():
    left = pd.DataFrame(
        {"chrom": ["chr1"], "start": [10], "end": [20], "name": ["left"]}
    )
    right = pd.DataFrame(
        {"chrom": ["chr1"], "start": [30], "end": [40], "name": ["right"]}
    )

    intersection = intersect_bins(left, right)

    assert intersection.empty
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


def test_intersect_bins_multiple_chromosomes():
    left = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [10, 10],
            "end": [30, 30],
            "name": ["left_1", "left_2"],
        }
    )
    right = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [20, 25],
            "end": [40, 35],
            "name": ["right_1", "right_2"],
        }
    )

    intersection = intersect_bins(left, right)

    assert intersection["chrom_1"].tolist() == ["chr1", "chr2"]
    assert intersection["name_2"].tolist() == ["right_1", "right_2"]
    assert intersection["overlap"].tolist() == [10, 5]


def test_intersect_bins_slack_expands_intervals():
    left = pd.DataFrame(
        {"chrom": ["chr1"], "start": [10], "end": [20], "name": ["left"]}
    )
    right = pd.DataFrame(
        {"chrom": ["chr1"], "start": [25], "end": [35], "name": ["right"]}
    )

    intersection = intersect_bins(left, right, slack=5)

    assert intersection.loc[0, "start_1"] == 5
    assert intersection.loc[0, "end_1"] == 25
    assert intersection.loc[0, "start_2"] == 20
    assert intersection.loc[0, "end_2"] == 40
    assert intersection.loc[0, "overlap"] == 5


@pytest.mark.parametrize(
    "strandedness,expected_names",
    [
        ("same", ["same"]),
        ("opposite", ["opposite"]),
    ],
)
def test_intersect_bins_strand_handling(strandedness, expected_names):
    left = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [10],
            "end": [30],
            "name": ["left"],
            "strand": ["+"],
        }
    )
    right = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [20, 20],
            "end": [40, 40],
            "name": ["same", "opposite"],
            "strand": ["+", "-"],
        }
    )

    intersection = intersect_bins(left, right, strandedness=strandedness)

    assert intersection["name_2"].tolist() == expected_names


def test_intersect_bins_synthesizes_unnamed_intervals():
    left = pd.DataFrame({"chrom": ["chr1"], "start": [10], "end": [30]})
    right = pd.DataFrame({"chrom": ["chr1"], "start": [20], "end": [40]})

    intersection = intersect_bins(left, right)

    assert intersection.loc[0, "name_1"] == "region_1_0"
    assert intersection.loc[0, "name_2"] == "region_2_0"
