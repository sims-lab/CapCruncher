import pytest
import pandas as pd
import pysam
import os

from ccanalyser.tools.filter import CCSliceFilter, TriCSliceFilter, TiledCSliceFilter

# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data")

def test_ccslice_filter():

    test_slices = os.path.join(dir_data, 'test', 'test_slices_to_filter_capture.tsv')
    df_test_slices = pd.read_csv(test_slices, sep=r'\s+')

    sf = CCSliceFilter(df_test_slices)
    n_removed = 0

    sf.remove_unmapped_slices()
    n_removed += 1 # 1 unmapped in test data
    assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)

    sf.remove_orphan_slices()
    n_removed += 1 # 1 orphan slice
    assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)

    sf.remove_blacklisted_slices()
    n_removed += 1 # 1 blacklisted slice
    assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)

    sf.remove_non_reporter_fragments()
    n_removed += 4 # 4 slices with no capture
    assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)








