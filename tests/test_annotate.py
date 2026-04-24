import os

import pandas as pd
import pytest
from capcruncher.api.annotate import annotate_intervals

# @pytest.fixture(scope="session")
# def ray_cluster():
#     import ray

#     ray.init(num_cpus=4, ignore_reinit_error=True, include_dashboard=False)


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "alignment_annotation")
    return data_dir


@pytest.mark.parametrize(
    "bed1,bed2,method,name,n_rows_expected",
    [
        (
            "test_slices_sorted.bed",
            "test_capture.bed",
            "count",
            "capture",
            4,
        ),
        (
            "test_slices_sorted.bed",
            "test_capture.bed",
            "get",
            "capture",
            4,
        ),
        (
            "test_slices_sorted.bed",
            "bad_bed.bed",
            "count",
            "capture_count",
            4,
        ),
        (
            "test_slices_sorted.bed",
            "blank.bed",
            "count",
            "capture_count",
            4,
        ),
    ],
)
def test_bed_intersection_succeeds(
    data_path, bed1, bed2, method, name, n_rows_expected
):

    intersection = annotate_intervals(
        query=os.path.join(data_path, bed1),
        annotations=os.path.join(data_path, bed2),
        name=name,
        method=method,
    )
    assert name in intersection.columns
    assert intersection.shape[0] == n_rows_expected


def test_bed_intersection_get_output(data_path):

    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    intersection = annotate_intervals(
        query=os.path.join(data_path, bed1),
        annotations=os.path.join(data_path, bed2),
        name="capture",
        method="get",
    )
    assert "capture" in intersection.columns
    assert intersection["capture"].value_counts().loc["CAPTURE"] == 1


def test_bed_intersection_count_output(data_path):
    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    intersection = annotate_intervals(
        query=os.path.join(data_path, bed1),
        annotations=os.path.join(data_path, bed2),
        name="capture",
        method="count",
    )
    assert "capture" in intersection.columns
    assert intersection["capture"].sum() == 1


def test_multi_intersection(data_path):
    slices = os.path.join(data_path, "test_slices_sorted.bed")
    capture = os.path.join(data_path, "test_capture.bed")
    capture_count = os.path.join(data_path, "test_capture_count.bed")
    blank = os.path.join(data_path, "blank.bed")

    for bed, action, name in zip(
        [capture, capture_count, blank],
        ["get", "count", "get"],
        ["capture", "capture_count", "blank"],
    ):
        slices = annotate_intervals(
            query=slices,
            annotations=bed,
            name=name,
            method=action,
        )


def test_get_intersection_preserves_input_row_order_and_count():
    slices = pd.DataFrame(
        {
            "chrom": ["chr2", "chr1", "chr1"],
            "start": [10, 100, 200],
            "end": [20, 200, 300],
            "name": ["first", "second", "third"],
        }
    )
    annotations = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [250, 150, 10],
            "end": [260, 175, 20],
            "name": ["ann_third", "ann_second", "ann_first"],
        }
    )

    intersection = annotate_intervals(
        query=slices,
        annotations=annotations,
        name="annotation",
        method="get",
    )

    assert intersection.shape[0] == slices.shape[0]
    assert intersection["Name"].tolist() == ["first", "second", "third"]
    assert intersection["annotation"].tolist() == [
        "ann_first",
        "ann_second",
        "ann_third",
    ]


def test_annotate_intervals_accepts_dataframe_inputs():
    slices = pd.DataFrame(
        {"chrom": ["chr1"], "start": [100], "end": [200], "name": ["slice"]}
    )
    annotations = pd.DataFrame(
        {"chrom": ["chr1"], "start": [100], "end": [150], "name": ["capture"]}
    )

    intersection = annotate_intervals(
        query=slices,
        annotations=annotations,
        name="annotation",
        method="get",
    )

    assert intersection["Name"].tolist() == ["slice"]
    assert intersection["annotation"].tolist() == ["capture"]


@pytest.mark.parametrize(
    "fraction, expected",
    [
        (0.5, "half"),
        (0.51, pd.NA),
    ],
)
def test_get_intersection_fraction_threshold_edges(fraction, expected):
    slices = pd.DataFrame(
        {"chrom": ["chr1"], "start": [100], "end": [200], "name": ["slice"]}
    )
    annotations = pd.DataFrame(
        {"chrom": ["chr1"], "start": [100], "end": [150], "name": ["half"]}
    )

    intersection = annotate_intervals(
        query=slices,
        annotations=annotations,
        name="annotation",
        method="get",
        fraction=fraction,
    )

    result = intersection["annotation"].iloc[0]
    if pd.isna(expected):
        assert pd.isna(result)
    else:
        assert result == expected


def test_get_intersection_numeric_annotation_names():
    slices = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [200, 300],
            "name": ["slice_a", "slice_b"],
        }
    )
    annotations = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [100],
            "end": [150],
            "name": [7],
        }
    )

    intersection = annotate_intervals(
        query=slices,
        annotations=annotations,
        name="annotation",
        method="get",
        fraction=0.5,
    )

    assert str(intersection["annotation"].dtype) == "Int64"
    assert intersection["annotation"].tolist() == [7, pd.NA]


def test_get_intersection_categorical_annotation_names():
    slices = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [100, 200],
            "end": [200, 300],
            "name": ["slice_a", "slice_b"],
        }
    )
    annotations = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [100],
            "end": [150],
            "name": pd.Categorical(["capture"]),
        }
    )

    intersection = annotate_intervals(
        query=slices,
        annotations=annotations,
        name="annotation",
        method="get",
        fraction=0.5,
    )

    assert isinstance(intersection["annotation"].dtype, pd.CategoricalDtype)
    assert intersection["annotation"].iloc[0] == "capture"
    assert pd.isna(intersection["annotation"].iloc[1])


@pytest.mark.parametrize("bed_file", ["blank.bed", "bad_bed.bed"])
def test_get_intersection_blank_and_bad_bed_fallback(data_path, bed_file):
    slices = os.path.join(data_path, "test_slices_sorted.bed")

    intersection = annotate_intervals(
        query=slices,
        annotations=os.path.join(data_path, bed_file),
        name="annotation",
        method="get",
    )

    assert intersection.shape[0] == 4
    assert "annotation" in intersection.columns
    assert intersection["annotation"].isna().all()
