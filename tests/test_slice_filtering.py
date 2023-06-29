import pytest
import os

from capcruncher.api.filter import CCSliceFilter, TriCSliceFilter, TiledCSliceFilter
from capcruncher.cli.alignments_filter import merge_annotations
from capcruncher.api.io import parse_bam


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "alignment_filtering")
    return data_dir


def get_slices(bam, annotations):
    df_alignment = parse_bam(bam)
    df_alignment = merge_annotations(df_alignment, annotations)
    return df_alignment


@pytest.mark.parametrize(
    "filter_class,bam,annotations,n_slices_expected",
    [
        (CCSliceFilter, "test.flashed.bam", "test.annotations.parquet", 135),
        (TriCSliceFilter, "test.flashed.bam", "test.annotations.parquet", 47),
        (TiledCSliceFilter, "test.flashed.bam", "test.annotations.parquet", 157),
    ],
)
def test_filters(data_path, filter_class, bam, annotations, n_slices_expected):

    bam = os.path.join(data_path, bam)
    annotations = os.path.join(data_path, annotations)

    df_slices = get_slices(bam, annotations)
    sf = filter_class(df_slices)

    sf.filter_slices()
    assert sf.slices.shape[0] == n_slices_expected
