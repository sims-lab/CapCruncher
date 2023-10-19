import pytest
import os
import pathlib
from typing import Union

from capcruncher.api.filter import CCSliceFilter, TriCSliceFilter, TiledCSliceFilter
from capcruncher.cli.alignments_filter import merge_annotations
from capcruncher.api.io import parse_bam


@pytest.fixture(scope="module")
def data_path():
    cwd = pathlib.Path.cwd()
    fn = pathlib.Path(__file__).relative_to(cwd)
    dirname = fn.parent
    data_dir = dirname / "data" / "alignment_filtering"
    return str(data_dir)


@pytest.fixture(scope="function")
def parquet_file(tmpdir):
    return pathlib.Path(tmpdir) / "test.parquet"


def get_slices(bam: str, annotations: str, parquet_file: Union[str, pathlib.Path]):
    df_alignment = parse_bam(bam)
    df_alignment.to_parquet(parquet_file)
    df_alignment = merge_annotations(parquet_file, annotations)
    return df_alignment


@pytest.mark.parametrize(
    "filter_class,bam,annotations,n_slices_expected",
    [
        (CCSliceFilter, "test.flashed.bam", "test.annotations.parquet", 135),
        (TriCSliceFilter, "test.flashed.bam", "test.annotations.parquet", 47),
        (TiledCSliceFilter, "test.flashed.bam", "test.annotations.parquet", 128),
    ],
)
def test_filters(
    data_path, filter_class, bam, annotations, n_slices_expected, parquet_file
):
    bam = os.path.join(data_path, bam)
    annotations = os.path.join(data_path, annotations)

    df_slices = get_slices(bam, annotations, parquet_file)
    sf = filter_class(df_slices)

    sf.filter_slices()
    assert sf.slices.shape[0] == n_slices_expected
