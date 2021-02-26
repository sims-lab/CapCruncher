import os

import pytest
import pandas as pd
from pybedtools import BedTool
from ccanalyser.tools.annotate import BedIntersection

# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data")


def test_bed_intersection():

    test_slices = os.path.join(dir_data, "test", "test_slices.bed")
    test_capture = os.path.join(dir_data, "test", "test_capture.bed")

    bt_slices = BedTool(test_slices)
    df_slices = bt_slices.to_dataframe()

    bi = BedIntersection(
        bed1=test_slices,
        bed2=test_capture,
        intersection_name="capture",
        intersection_method="count",
        n_cores=1,
    )

    intersection = bi.intersection
    assert intersection.name == 'capture'
    assert intersection.shape[0] == df_slices.shape[0]