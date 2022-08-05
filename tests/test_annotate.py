import os

import numpy as np
import pandas as pd
import pytest
from capcruncher.tools.annotate import BedIntersection
from more_itertools import intersperse
from pandas.api.types import is_numeric_dtype, is_string_dtype
from pybedtools import BedTool


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "alignment_annotation")
    return data_dir


@pytest.mark.parametrize(
    "bed1,bed2,method,n_rows_expected,dtype,dtype_check_func",
    [
        (
            "test_slices_sorted.bed",
            "test_capture.bed",
            "count",
            4,
            "Int8",
            is_numeric_dtype,
        ),
        (
            "test_slices_sorted.bed",
            "test_capture.bed",
            "get",
            1,
            "string",
            is_string_dtype,
        ),
        ("test_slices_sorted.bed", "bad_bed.bed", "count", 4, "Int8", is_numeric_dtype),
        ("test_slices_sorted.bed", "blank.bed", "count", 4, "Int8", is_numeric_dtype),
    ],
)
def test_bed_intersection_succeeds(
    data_path, bed1, bed2, method, n_rows_expected, dtype, dtype_check_func
):

    bi = BedIntersection(
        bed1=os.path.join(data_path, bed1),
        bed2=os.path.join(data_path, bed2),
        intersection_name="capture",
        intersection_method=method,
        invalid_bed_action="ignore",
        dtype=dtype,
    )

    intersection = bi.get_intersection()
    assert intersection.name == "capture"
    assert intersection.shape[0] == n_rows_expected
    assert dtype_check_func(intersection.values)


@pytest.mark.parametrize(
    "bed1,bed2,method",
    [
        (
            "test_slices.bed",
            "test_capture.bed",
            "count",
        ),
        (
            "test_slices_sorted.bed",
            "bad_bed.bed",
            "count",
        ),
        (
            "test_slices_sorted.bed",
            "blank.bed",
            "count",
        ),
    ],
)
@pytest.mark.xfail
def test_bed_intersection_fails(data_path, bed1, bed2, method):

    bi = BedIntersection(
        bed1=os.path.join(data_path, bed1),
        bed2=os.path.join(data_path, bed2),
        intersection_name="capture",
        intersection_method=method,
        invalid_bed_action="error",
    )

    intersection = bi.get_intersection()


def test_bed_intersection_count_succeeds(data_path):
    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    bi = BedIntersection(
        bed1=os.path.join(data_path, bed1),
        bed2=os.path.join(data_path, bed2),
        intersection_name="capture",
        intersection_method="count",
        invalid_bed_action="error",
        dtype="Int8",
    )

    intersection = bi.get_intersection()
    assert intersection.name == "capture"
    assert intersection.shape[0] == 4
    assert is_numeric_dtype(intersection.values)
    assert intersection.sum() == 1


def test_bed_intersection_get_succeeds(data_path):
    bed1 = "test_slices_sorted.bed"
    bed2 = "test_capture.bed"

    bi = BedIntersection(
        bed1=os.path.join(data_path, bed1),
        bed2=os.path.join(data_path, bed2),
        intersection_name="capture",
        intersection_method="get",
        invalid_bed_action="error",
        dtype="category",
    )

    intersection = bi.get_intersection()
    assert intersection.name == "capture"
    assert intersection.value_counts().loc["CAPTURE"] == 1
