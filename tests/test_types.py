import json
import os
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from capcruncher.api.alignments.annotate import AlignmentAnnotateOptions, annotate
from capcruncher.api.alignments.filter import AlignmentFilterOptions
from capcruncher.api.interactions.compare import concat, summarise
from capcruncher.cli import cli
from capcruncher.types import AnnotationAction, Assay, DuplicateAction, ReadType
from capcruncher.utils import load_dict, save_dict

_DATA_DIR = Path(os.path.dirname(__file__)) / "data" / "alignment_annotation"


def test_alignment_annotate_options_use_string_enums(tmp_path):
    slices = tmp_path / "slices.bed"
    targets = tmp_path / "targets.bed"
    slices.touch()
    targets.touch()
    options = AlignmentAnnotateOptions(
        slices=slices,
        actions=["get"],
        bed_files=[targets],
        names=["targets"],
        duplicates="remove",
    )

    assert options.actions == (AnnotationAction.GET,)
    assert options.duplicates is DuplicateAction.REMOVE

    with pytest.raises(ValueError, match="actions must be one of: get, count"):
        AlignmentAnnotateOptions(
            slices=slices,
            actions=["fetch"],
            bed_files=[targets],
            names=["targets"],
        )

    with pytest.raises(ValueError, match="duplicates must be one of: remove"):
        AlignmentAnnotateOptions(slices=slices, duplicates="keep")


def test_annotate_accepts_path_input_for_bam(tmp_path):
    bam_path = _DATA_DIR / "test.pe.bam"
    output = tmp_path / "annotated.parquet"

    annotate(slices=bam_path, output=output)

    assert output.exists()


def test_alignment_filter_options_validate_string_enums(tmp_path):
    bam = tmp_path / "reads.bam"
    annotations = tmp_path / "annotations.parquet"
    bam.touch()
    annotations.touch()

    options = AlignmentFilterOptions(
        bam=bam,
        annotations=annotations,
        method="capture",
        read_type="flashed",
    )

    assert options.method is Assay.CAPTURE
    assert options.read_type is ReadType.FLASHED

    with pytest.raises(ValueError, match="method must be one of: capture, tri, tiled"):
        AlignmentFilterOptions(bam=bam, annotations=annotations, method="invalid")

    with pytest.raises(ValueError, match="read_type must be one of: flashed, pe"):
        AlignmentFilterOptions(bam=bam, annotations=annotations, read_type="PE")


def test_alignment_filter_cli_rejects_invalid_string_enum():
    result = CliRunner().invoke(
        cli,
        [
            "alignments",
            "filter",
            "invalid",
            "--bam",
            "reads.bam",
            "--annotations",
            "annotations.parquet",
        ],
    )

    assert result.exit_code != 0
    assert "invalid" in result.output


def test_summarise_rejects_unsupported_method(tmp_path):
    infile = tmp_path / "union.tsv"
    infile.write_text("chrom\tstart\tend\tsample\nchr1\t1\t2\t3\n")

    with pytest.raises(ValueError, match="summary_methods must be one of: mean"):
        summarise(infile=infile, summary_methods=("median",))


def test_compare_concat_and_summarise_use_polars_outputs(tmp_path):
    bedgraph_a = tmp_path / "a.bedgraph"
    bedgraph_b = tmp_path / "b.bedgraph"
    union_path = tmp_path / "union.tsv"

    bedgraph_a.write_text("chr1\t0\t10\t2\nchr1\t0\t10\t3\nchr1\t10\t20\t4\n")
    bedgraph_b.write_text("chr1\t0\t10\t7\nchr1\t20\t30\t5\n")

    union_by_viewpoint = concat(
        [bedgraph_a, bedgraph_b],
        viewpoint="vp1",
        format="bedgraph",
        output=union_path,
    )

    union = union_by_viewpoint["vp1"]
    assert isinstance(union, pl.DataFrame)
    assert union.to_dicts() == [
        {"chrom": "chr1", "start": 0, "end": 10, "a.bedgraph": 5, "b.bedgraph": 7},
        {"chrom": "chr1", "start": 10, "end": 20, "a.bedgraph": 4, "b.bedgraph": 0},
        {"chrom": "chr1", "start": 20, "end": 30, "a.bedgraph": 0, "b.bedgraph": 5},
    ]

    summarise(
        infile=union_path,
        output_prefix=tmp_path / "summary.",
        group_names=("grp",),
        group_columns=("0,1",),
    )

    summary = pl.read_csv(
        tmp_path / "summary.grp.mean-summary.bedgraph",
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "score"],
    )
    assert summary.to_dicts() == [
        {"chrom": "chr1", "start": 0, "end": 10, "score": 6.0},
        {"chrom": "chr1", "start": 10, "end": 20, "score": 2.0},
        {"chrom": "chr1", "start": 20, "end": 30, "score": 2.5},
    ]


def test_dict_serialisation_rejects_invalid_format_and_dtype(tmp_path):
    json_path = tmp_path / "mapping.json"
    save_dict({1: 2}, json_path, format="json")
    assert json.loads(json_path.read_text()) == {"1": 2}

    with pytest.raises(ValueError, match="format must be one of: json, pickle"):
        save_dict({1: 2}, tmp_path / "mapping.bad", format="toml")

    with pytest.raises(ValueError, match="dtype must be one of: int, str"):
        load_dict(json_path, format="json", dtype="float")
