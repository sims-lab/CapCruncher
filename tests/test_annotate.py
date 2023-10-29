import os

import pytest
from capcruncher.api.annotate import BedIntersector

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

    bi = BedIntersector(
        bed_a=os.path.join(data_path, bed1),
        bed_b=os.path.join(data_path, bed2),
        name=name,
    )

    intersection = bi.get_intersection(method=method)
    assert name in intersection.columns
    assert intersection.df.shape[0] == n_rows_expected


def test_bed_intersection_get_output(data_path):

    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    bi = BedIntersector(
        bed_a=os.path.join(data_path, bed1),
        bed_b=os.path.join(data_path, bed2),
        name="capture",
    )

    intersection = bi.get_intersection(method="get")
    assert "capture" in intersection.columns
    assert intersection.df["capture"].value_counts().loc["CAPTURE"] == 1


def test_bed_intersection_count_output(data_path):

    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    bi = BedIntersector(
        bed_a=os.path.join(data_path, bed1),
        bed_b=os.path.join(data_path, bed2),
        name="capture",
    )

    intersection = bi.get_intersection(method="count")
    assert "capture" in intersection.columns
    assert intersection.df["capture"].sum() == 1
