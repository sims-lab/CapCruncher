# import pytest
# import pandas as pd
# import pysam
# import os
# import numpy as np

# from capcruncher.tools.filter import CCSliceFilter, TriCSliceFilter, TiledCSliceFilter

# # Pre-run setup
# dir_test = os.path.realpath(os.path.dirname(__file__))
# dir_package = os.path.dirname(dir_test)
# dir_data = os.path.join(dir_package, "data")

# test_slices = os.path.join(dir_data, 'test', 'test_slices_to_filter_capture.tsv')
# df_test_slices = pd.read_csv(test_slices, sep='\t')
# df_test_slices['restriction_fragment'] = df_test_slices['restriction_fragment'].replace('.', -1).astype(int)

# test_yaml = os.path.join(dir_data, 'test', 'ccslicefilter_test.yml')

# def test_ccslicefilter():

#     sf = CCSliceFilter(df_test_slices)
#     n_removed = 0

#     sf.remove_unmapped_slices()
#     n_removed += 1 # 1 unmapped in test data
#     assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)

#     sf.remove_orphan_slices()
#     n_removed += 1 # 1 orphan slice
#     assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)

#     sf.remove_blacklisted_slices()
#     n_removed += 1 # 1 blacklisted slice
#     assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)

#     sf.remove_non_reporter_fragments()
#     n_removed += 4 # 4 slices with no capture
#     assert sf.slices.shape[0] == (df_test_slices.shape[0] - n_removed)

# def test_ccslicefilter_filter_slices():

#     sf = CCSliceFilter(df_test_slices)

#     sf.filter_slices()
#     assert sf.slices.shape[0] == 2

# def test_from_yaml():

#     sf = CCSliceFilter(df_test_slices, filter_stages=test_yaml)

#     sf.filter_slices()
#     assert sf.slices.shape[0] == 2







