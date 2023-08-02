import os

import pytest
from capcruncher.api.annotate import BedFileIntersection
from pandas.api.types import (
    is_numeric_dtype,
    is_categorical_dtype,
    is_object_dtype,
)
import ray


@pytest.fixture(scope="session")
def ray_cluster():
    import ray

    ray.init(num_cpus=4, ignore_reinit_error=True, include_dashboard=False)


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "alignment_annotation")
    return data_dir


@pytest.mark.parametrize(
    "bed1,bed2,method,name,n_rows_expected,dtype_check_func",
    [
        (
            "test_slices_sorted.bed",
            "test_capture.bed",
            "count",
            "capture",
            1,
            is_numeric_dtype,
        ),
        (
            "test_slices_sorted.bed",
            "test_capture.bed",
            "get",
            "capture",
            1,
            is_categorical_dtype,
        ),
        (
            "test_slices_sorted.bed",
            "bad_bed.bed",
            "count",
            "capture_count",
            4,
            is_object_dtype,
        ),
        (
            "test_slices_sorted.bed",
            "blank.bed",
            "count",
            "capture_count",
            4,
            is_object_dtype,
        ),
    ],
)
def test_bed_intersection_succeeds(
    ray_cluster, data_path, bed1, bed2, method, name, n_rows_expected, dtype_check_func
):

    bi = BedFileIntersection.remote(
        bed_a=os.path.join(data_path, bed1),
        bed_b=os.path.join(data_path, bed2),
        name=name,
        action=method,
    )

    intersection = ray.get(bi.intersection.remote())
    assert intersection.name == name
    assert intersection.shape[0] == n_rows_expected
    assert dtype_check_func(intersection.values)


def test_bed_intersection_get_output(data_path):

    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    bi = BedFileIntersection.remote(
        bed_a=os.path.join(data_path, bed1),
        bed_b=os.path.join(data_path, bed2),
        name="capture",
        action="get",
    )

    intersection = ray.get(bi.intersection.remote())
    assert intersection.name == "capture"
    assert intersection.value_counts().loc["CAPTURE"] == 1


def test_bed_intersection_count_output(data_path):

    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    bi = BedFileIntersection.remote(
        bed_a=os.path.join(data_path, bed1),
        bed_b=os.path.join(data_path, bed2),
        name="capture",
        action="count",
    )

    intersection = ray.get(bi.intersection.remote())
    assert intersection.name == "capture"
    assert intersection.sum() == 1
