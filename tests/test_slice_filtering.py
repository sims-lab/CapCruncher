import pytest
import os
import pathlib

import pandas as pd
import polars as pl

from capcruncher.api.filtering.pipeline import CCSliceFilter, TriCSliceFilter, TiledCSliceFilter
from capcruncher.api.alignments.filter import merge_annotations, remove_unused_categories
from capcruncher.api.alignments.io import parse_bam


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


def get_slices(bam: str, annotations: str, parquet_file: str | pathlib.Path):
    df_alignment = parse_bam(bam)
    assert isinstance(df_alignment, pl.DataFrame)
    df_alignment.write_parquet(parquet_file)
    df_alignment = merge_annotations(parquet_file, annotations)
    return df_alignment


def test_merge_annotations_normalises_join_key_dtypes(tmp_path):
    slices = tmp_path / "slices.parquet"
    annotations = tmp_path / "annotations.parquet"

    pl.DataFrame(
        {
            "slice_name": ["slice-a"],
            "chrom": ["chr1"],
            "start": [10],
            "slice_id": [1],
        },
        schema_overrides={"chrom": pl.Categorical},
    ).write_parquet(slices)
    pl.DataFrame(
        {
            "slice_name": ["slice-a"],
            "Chromosome": ["chr1"],
            "Start": [10],
            "End": [20],
            "capture": ["vp1"],
        }
    ).write_parquet(annotations)

    df_alignment = merge_annotations(slices, annotations)

    assert isinstance(df_alignment, pl.DataFrame)
    assert df_alignment[["slice_name", "chrom", "start", "capture"]].to_dicts() == [
        {
            "slice_name": "slice-a",
            "chrom": "chr1",
            "start": 10,
            "capture": "vp1",
        }
    ]


def test_remove_unused_categories_prunes_viewpoint_labels():
    df = pd.DataFrame(
        {
            "viewpoint": pd.Categorical(
                ["Slc25A37", "Slc25A37"],
                categories=["Slc25A37", "reporters_pe_80", "duplicate_coords_1"],
            )
        }
    )

    cleaned = remove_unused_categories(df)

    assert cleaned["viewpoint"].cat.categories.to_list() == ["Slc25A37"]


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


def test_default_filter_stage_stats(data_path, parquet_file):
    bam = os.path.join(data_path, "test.flashed.bam")
    annotations = os.path.join(data_path, "test.annotations.parquet")
    df_slices = get_slices(bam, annotations, parquet_file)

    sf = CCSliceFilter(df_slices)
    sf.filter_slices()

    assert [
        (stat.stage, stat.n_fragments, stat.n_slices) for stat in sf.filtering_stats
    ] == [
        ("pre-filtering", 91, 192),
        ("mapped", 91, 192),
        ("contains_single_capture", 78, 179),
        ("contains_capture_and_reporter", 68, 157),
        ("duplicate_filtered", 59, 135),
    ]


def test_toml_filter_profile_controls_stage_order(data_path, parquet_file, tmp_path):
    profile = tmp_path / "filter_profile.toml"
    profile.write_text(
        """
assay = "capture"

[[stages]]
name = "pre-filtering"
steps = ["get_unfiltered_slices"]

[[stages]]
name = "mapped"
steps = ["remove_unmapped_slices"]
""".strip()
    )
    bam = os.path.join(data_path, "test.flashed.bam")
    annotations = os.path.join(data_path, "test.annotations.parquet")
    df_slices = get_slices(bam, annotations, parquet_file)

    sf = CCSliceFilter(df_slices, filter_profile=profile)
    sf.filter_slices()

    assert sf.slices.shape[0] == 192
    assert [stat.stage for stat in sf.filtering_stats] == ["pre-filtering", "mapped"]


def test_toml_filter_profile_rejects_unknown_step(tmp_path):
    profile = tmp_path / "filter_profile.toml"
    profile.write_text(
        """
assay = "capture"

[[stages]]
name = "bad"
steps = ["remove_everything"]
""".strip()
    )

    with pytest.raises(ValueError, match="Unknown filter step"):
        CCSliceFilter(pd.DataFrame(), filter_profile=profile)


def test_toml_filter_profile_rejects_duplicate_stage_names(tmp_path):
    profile = tmp_path / "filter_profile.toml"
    profile.write_text(
        """
assay = "capture"

[[stages]]
name = "mapped"
steps = ["get_unfiltered_slices"]

[[stages]]
name = "mapped"
steps = ["remove_unmapped_slices"]
""".strip()
    )

    with pytest.raises(ValueError, match="Duplicate filter stage name"):
        CCSliceFilter(pd.DataFrame(), filter_profile=profile)


def test_toml_filter_profile_rejects_wrong_assay(tmp_path):
    profile = tmp_path / "filter_profile.toml"
    profile.write_text(
        """
assay = "tiled"

[[stages]]
name = "pre-filtering"
steps = ["get_unfiltered_slices"]
""".strip()
    )

    with pytest.raises(ValueError, match="does not match requested assay"):
        CCSliceFilter(pd.DataFrame(), filter_profile=profile)


def test_toml_filter_profile_rejects_assay_incompatible_step(tmp_path):
    profile = tmp_path / "filter_profile.toml"
    profile.write_text(
        """
assay = "capture"

[[stages]]
name = "bad"
steps = ["remove_religation"]
""".strip()
    )

    with pytest.raises(ValueError, match="is not valid for assay"):
        CCSliceFilter(pd.DataFrame(), filter_profile=profile)


def test_yaml_filter_profiles_are_not_supported(tmp_path):
    profile = tmp_path / "filter_profile.yml"
    profile.write_text("pre-filtering:\n  - get_unfiltered_slices\n")

    with pytest.raises(ValueError, match="YAML filter profiles are no longer supported"):
        CCSliceFilter(pd.DataFrame(), filter_profile=profile)
