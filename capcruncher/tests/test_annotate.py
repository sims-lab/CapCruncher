import os

import pytest
import pandas as pd
from pybedtools import BedTool
from capcruncher.tools.annotate import BedIntersection
import numpy as np
from pandas.api.types import is_numeric_dtype, is_string_dtype

# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data", "test", "alignment_annotation")


@pytest.mark.parametrize(
    "bed1,bed2,method,n_rows_expected,dtype_func",
    [
        (
            os.path.join(dir_data, "test_slices_sorted.bed"),
            os.path.join(dir_data, "test_capture.bed"),
            "count",
            4,
            is_numeric_dtype
        ),
        (
            os.path.join(dir_data, "test_slices_sorted.bed"),
            os.path.join(dir_data, "test_capture.bed"),
            "get",
            4,
            is_string_dtype
        ),
        (
            os.path.join(dir_data, "test_slices_sorted.bed"),
            os.path.join(dir_data, "bad_bed.bed"),
            "count",
            4,
            is_numeric_dtype
        ),
        (
            os.path.join(dir_data, "test_slices_sorted.bed"),
            os.path.join(dir_data, "blank.bed"),
            "count",
            4,
            is_numeric_dtype
        ),
    ],
)
def test_bed_intersection_succeeds(bed1, bed2, method, n_rows_expected, dtype_func):

    bi = BedIntersection(
        bed1=bed1,
        bed2=bed2,
        intersection_name="capture",
        intersection_method=method,
        invalid_bed_action="ignore",
    )

    intersection = bi.get_intersection()
    assert intersection.name == "capture"
    assert intersection.shape[0] == n_rows_expected
    assert dtype_func(intersection.values)




@pytest.mark.parametrize(
    "bed1,bed2,method",
    [
        (
            os.path.join(dir_data, "test_slices.bed"),
            os.path.join(dir_data, "test_capture.bed"),
            "count",
        ),
        (
            os.path.join(dir_data, "test_slices_sorted.bed"),
            os.path.join(dir_data, "bad_bed.bed"),
            "count",
        ),
        (
            os.path.join(dir_data, "test_slices_sorted.bed"),
            os.path.join(dir_data, "blank.bed"),
            "count",
        ),
    ],
)
@pytest.mark.xfail
def test_bed_intersection_fails(bed1, bed2, method):

    bi = BedIntersection(
        bed1=bed1,
        bed2=bed2,
        intersection_name="capture",
        intersection_method=method,
        invalid_bed_action="error",
    )

    intersection = bi.get_intersection()






# def test_bed_intersection_bad_format():

#     test_slices = os.path.join(dir_data, "test", "test_slices.bed")
#     test_capture = os.path.join(dir_data, "test", "test_capture_bad_format.bed")

#     bt_slices = BedTool(test_slices)
#     df_slices = bt_slices.to_dataframe()

#     #Without ignoring error
#     try:

#         bi = BedIntersection(
#             bed1=test_slices,
#             bed2=test_capture,
#             intersection_name="capture",
#             intersection_method="count",
#             )

#         intersection = bi.intersection

#     except Exception as e:
#         assert isinstance(e, ValueError)


#     #With ignoring error

#     bi = BedIntersection(
#         bed1=test_slices,
#         bed2=test_capture,
#         intersection_name="capture",
#         intersection_method="count",
#         invalid_bed_action="ignore")

#     intersection = bi.intersection
#     assert intersection.name == 'capture'
#     assert intersection.shape[0] == df_slices.shape[0]
#     assert is_numeric_dtype(intersection.values)


# def test_bed_intersection_blank():

#     test_slices = os.path.join(dir_data, "test", "test_slices.bed")
#     test_capture = os.path.join(dir_data, "test", "test_capture_blank.bed")

#     bt_slices = BedTool(test_slices)
#     df_slices = bt_slices.to_dataframe()

#     bi = BedIntersection(
#         bed1=test_slices,
#         bed2=test_capture,
#         intersection_name="capture",
#         intersection_method="count",
#         invalid_bed_action="ignore")

#     intersection = bi.intersection
#     assert intersection.name == 'capture'
#     assert intersection.shape[0] == df_slices.shape[0]
#     assert is_numeric_dtype(intersection.values)
