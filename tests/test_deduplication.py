import glob
import os
import pickle

import pandas as pd
import pytest
import ujson
from capcruncher.cli import fastq_deduplicate
from capcruncher.tools.deduplicate import (
    ReadDeduplicationParserProcess,
    ReadDuplicateRemovalProcess,
)
from capcruncher.tools.io import FastqReaderProcess
from capcruncher.utils import MockFastqRecord, get_file_type, load_dict, save_dict
from pybedtools.bedtool import BedTool
import pysam


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "fastq_deduplication")
    return data_dir


@pytest.mark.parametrize(
    "fastq_files,outfile,format,n_records_expected",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "parsed.json",
            "json",
            1520,
        ),
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "parsed.pkl",
            "pickle",
            1520,
        ),
    ],
)
def test_fastq_parsing(
    data_path, tmpdir, fastq_files, outfile, format, n_records_expected
):

    infiles = [os.path.join(data_path, fn) for fn in fastq_files]

    outfile_path = os.path.join(tmpdir, outfile)

    fastq_deduplicate.parse(infiles, outfile_path)

    assert os.path.exists(outfile_path)

    with open(outfile_path, "rb") as output:
        if format == "json":
            parse_test = ujson.load(output)
        elif format == "pickle":
            parse_test = pickle.load(output)

    assert len(parse_test) == n_records_expected


@pytest.mark.parametrize(
    "infiles,outfile,n_duplicates_expected",
    [
        (("parsed.json",), "deduplicated.json", 538),
        (("parsed.pickle",), "deduplicated.pickle", 538),
        (("parsed.json", "parsed.pickle"), "deduplicated.pickle", 1520),
    ],
)
def test_fastq_duplicate_identification(
    data_path, tmpdir, infiles, outfile, n_duplicates_expected
):
    infiles_paths = [os.path.join(data_path, fn) for fn in infiles]
    outfile_path = os.path.join(tmpdir, outfile)
    duplicated_ids = fastq_deduplicate.identify(infiles_paths, outfile_path)

    assert os.path.exists(outfile_path)

    outfile_type = get_file_type(outfile)
    output_dict = load_dict(outfile_path, format=outfile_type)

    assert len(duplicated_ids) == n_duplicates_expected


@pytest.mark.parametrize(
    "infiles,duplicates,prefix,n_duplicates_expected",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "identified.pkl",
            "out",
            538,
        )
    ],
)
def test_fastq_duplicate_removal(
    data_path, tmpdir, infiles, duplicates, prefix, n_duplicates_expected
):

    infiles_paths = [os.path.join(data_path, fn) for fn in infiles]
    ids = os.path.join(data_path, duplicates)
    out_prefix = os.path.join(tmpdir, prefix)

    stats = fastq_deduplicate.remove(
        infiles_paths,
        duplicated_ids=ids,
        output_prefix=out_prefix,
        stats_prefix=os.path.join(tmpdir, "dedup"),
    )

    stats = stats.set_index("stat_type")

    assert stats.at["reads_removed", "stat"] == n_duplicates_expected

    outfiles = glob.glob(f"{out_prefix}*")
    assert len(outfiles) == len(infiles)

    assert (
        len([r for r in pysam.FastxFile(outfiles[0])]) == stats.at["reads_unique", "stat"])

