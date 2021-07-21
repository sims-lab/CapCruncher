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
dir_data = os.path.join(dir_package, "data")

def test_bed_intersection_correct():

    test_slices = os.path.join(dir_data, "test", "test_slices.bed")
    test_capture = os.path.join(dir_data, "test", "test_capture.bed")

    bt_slices = BedTool(test_slices)
    df_slices = bt_slices.to_dataframe()

    # Test "count" method
    bi = BedIntersection(
        bed1=test_slices,
        bed2=test_capture,
        intersection_name="capture",
        intersection_method="count",
    )

    intersection = bi.intersection
    assert intersection.name == 'capture'
    assert intersection.shape[0] == df_slices.shape[0]
    assert is_numeric_dtype(intersection.values)
    

    # Test "get" method
    bi = BedIntersection(
        bed1=test_slices,
        bed2=test_capture,
        intersection_name="capture",
        intersection_method="get",
    )

    intersection = bi.intersection
    assert intersection.name == 'capture'
    assert intersection.shape[0] == df_slices.shape[0]
    assert is_string_dtype(intersection.values)

def test_bed_intersection_bad_format():

    test_slices = os.path.join(dir_data, "test", "test_slices.bed")
    test_capture = os.path.join(dir_data, "test", "test_capture_bad_format.bed")

    bt_slices = BedTool(test_slices)
    df_slices = bt_slices.to_dataframe()

    #Without ignoring error
    try:

        bi = BedIntersection(
            bed1=test_slices,
            bed2=test_capture,
            intersection_name="capture",
            intersection_method="count",
            n_cores=1)
        
        intersection = bi.intersection

    except Exception as e:
        assert isinstance(e, ValueError)
    

    #With ignoring error

    bi = BedIntersection(
        bed1=test_slices,
        bed2=test_capture,
        intersection_name="capture",
        intersection_method="count",
        n_cores=1,
        invalid_bed_action="ignore")
    
    intersection = bi.intersection
    assert intersection.name == 'capture'
    assert intersection.shape[0] == df_slices.shape[0]
    assert is_numeric_dtype(intersection.values)


def test_bed_intersection_blank():

    test_slices = os.path.join(dir_data, "test", "test_slices.bed")
    test_capture = os.path.join(dir_data, "test", "test_capture_blank.bed")

    bt_slices = BedTool(test_slices)
    df_slices = bt_slices.to_dataframe()
    
    bi = BedIntersection(
        bed1=test_slices,
        bed2=test_capture,
        intersection_name="capture",
        intersection_method="count",
        n_cores=1,
        invalid_bed_action="ignore")
    
    intersection = bi.intersection
    assert intersection.name == 'capture'
    assert intersection.shape[0] == df_slices.shape[0]
    assert is_numeric_dtype(intersection.values)






    