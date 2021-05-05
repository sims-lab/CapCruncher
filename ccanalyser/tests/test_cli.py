import glob
import os

import click
import pysam
import pytest
import ujson
import xopen
from click.testing import CliRunner
from ccanalyser.cli import cli
import pandas as pd


# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data")


@pytest.fixture(scope="session", autouse=True)
def cleanup():

    yield

    # Remove previous test files
    for fn in glob.glob(os.path.join(dir_test, "test", "*")):
        os.unlink(fn)

    for fn in glob.glob(os.path.join(dir_test, "stats", "*")):
        os.unlink(fn)


def get_fastq_n_records(fq):
    return sum(1 for l in pysam.FastxFile(fq))


class DataDeduplicate:
    fq1 = os.path.join(dir_data, "test", "duplicated_1.fastq.gz")
    fq2 = os.path.join(dir_data, "test", "duplicated_2.fastq.gz")
    result_parsed = os.path.join(dir_test, "expected", "fq_parsed_result.json")
    result_identify = os.path.join(dir_test, "expected", "fq_duplicates_results.json")
    result_remove = os.path.join(dir_test, "expected", "fq_dedup_1.fastq")


# Test data
data_dd = DataDeduplicate


def test_cli_runs():
    """Test checks that the cli is functional and the help option works"""

    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0


def test_genome_digest():

    test_fa = os.path.join(dir_data, "test", "test.fa")
    test_output = os.path.join(dir_test, "test", "test_digest_genome.bed")
    test_output_stats = os.path.join(dir_test, "stats", "test_digest_genome.bed")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "genome",
            "digest",
            test_fa,
            "-r",
            "dpnii",
            "-o",
            test_output,
            "-l",
            test_output_stats,
            "--sort",
        ],
    )

    assert result.exit_code == 0

    with open(test_output) as r:
        test_bed = [l.strip().split("\t") for l in r]

    # Should be just the one cutsite at the moment
    assert len(test_bed) == 2


def test_deduplicate_parse():

    # Test parsing
    output_parsed = os.path.join(dir_test, "test", "fq_parsed_test.json")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "fastq",
            "deduplicate",
            "parse",
            data_dd.fq1,
            data_dd.fq2,
            "-o",
            output_parsed,
        ],
    )

    # Check that the script exits successfully
    assert result.exit_code == 0

    with open(output_parsed) as f:
        result_test = ujson.load(f)

    with open(data_dd.result_parsed) as f:
        result_correct = ujson.load(f)

    # Check that all keys and values match between the saved file and the test
    assert result_test == result_correct


def test_fastq_deduplicate_identification():
    output_duplicates_test = os.path.join(dir_test, "test", "fq_duplicates_test.json")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "fastq",
            "deduplicate",
            "identify",
            data_dd.result_parsed,
            "-o",
            output_duplicates_test,
        ],
    )

    assert result.exit_code == 0

    with open(output_duplicates_test) as f:
        result_test = ujson.load(f)

    with open(data_dd.result_identify) as f:
        result_correct = ujson.load(f)

    # Checks that the same number of duplicates are identified, some randomness in identification
    assert len(result_test) == len(result_correct)


def test_fastq_deduplicate_removal():

    # Test removal
    output_removal_prefix = os.path.join(dir_test, "test", "fq_dedup_test")
    output_removal_test = output_removal_prefix + "_1.fastq"
    output_removal_test_stats = os.path.join(dir_test, "stats", "deduplication")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "fastq",
            "deduplicate",
            "remove",
            "-d",
            data_dd.result_identify,
            "-o",
            output_removal_prefix,
            data_dd.fq1,
            data_dd.fq2,
            "--stats_prefix",
            output_removal_test_stats,
            "--sample_name",
            "test",
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(output_removal_test_stats + ".deduplication.csv")

    with open(output_removal_test) as r:
        result_test = r.readlines()

    with open(data_dd.result_remove) as r:
        result_correct = r.readlines()

    with open(data_dd.result_identify) as r:
        duplicates = ujson.load(r)

    fq_unfilt_n_entries = get_fastq_n_records(data_dd.fq1)
    fq_dd_n_entries = get_fastq_n_records(output_removal_test)

    # Checks the number of expected duplicates are removed
    assert len(duplicates) == fq_unfilt_n_entries - fq_dd_n_entries

    # Checks the test file matches the expected file
    assert len(result_test) == len(result_correct)


def test_fastq_digest():

    fq1 = os.path.join(dir_data, "test", "digest_1.fastq.gz")
    fq2 = os.path.join(dir_data, "test", "digest_2.fastq.gz")
    test_output_flashed = os.path.join(
        dir_test, "test", "test_fastq_digest_flashed.fastq"
    )
    test_output_pe = os.path.join(dir_test, "test", "test_fastq_digest_pe.fastq")
    test_output_stats = os.path.join(dir_test, "stats", "digestion")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "fastq",
            "digest",
            "-m",
            "pe",
            "-r",
            "dpnii",
            "-o",
            test_output_pe,
            "--stats_prefix",
            test_output_stats,
            "--sample_name",
            "test",
            fq1,
            fq2,
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(test_output_stats + ".digestion.read.summary.csv")

    os.unlink(test_output_stats + ".digestion.read.summary.csv")
    result = runner.invoke(
        cli,
        [
            "fastq",
            "digest",
            "-m",
            "flashed",
            "-r",
            "dpnii",
            "-o",
            test_output_flashed,
            "--stats_prefix",
            test_output_stats,
            "--sample_name",
            "test",
            fq1,
        ],
    )

    assert result.exit_code == 0

    # check stats are produced
    assert os.path.exists(test_output_stats + ".digestion.read.summary.csv")


def test_alignments_annotate():

    test_output = os.path.join(dir_test, "test", "test_annotate.tsv")
    test_bed = os.path.join(dir_data, "test", "test_slices.bed")
    test_capture = os.path.join(dir_data, "test", "test_capture.bed")
    test_exclusion = os.path.join(dir_data, "test", "test_exclusions.bed")
    test_rf = os.path.join(dir_data, "test", "test_rf.bed")

    # Test with bad exclusions file
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "alignments",
            "annotate",
            test_bed,
            "-b",
            test_capture,
            "-b",
            test_capture,
            "-b",
            test_exclusion,
            "-b",
            test_rf,
            "-a",
            "get",
            "-a",
            "count",
            "-a",
            "get",
            "-a",
            "get",
            "-n",
            "capture",
            "-n",
            "capture_count",
            "-n",
            "exclusion",
            "-n",
            "rf",
            "-o",
            test_output,
        ],
    )

    assert result.exit_code == 1

    # Test with ignoring bad file
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "alignments",
            "annotate",
            test_bed,
            "-b",
            test_capture,
            "-b",
            test_capture,
            "-b",
            test_exclusion,
            "-b",
            test_rf,
            "-a",
            "get",
            "-a",
            "count",
            "-a",
            "get",
            "-a",
            "get",
            "-n",
            "capture",
            "-n",
            "capture_count",
            "-n",
            "exclusion",
            "-n",
            "rf",
            "-o",
            test_output,
            "--invalid_bed_action",
            "ignore",
        ],
    )

    assert result.exit_code == 0

    # Test with bad path
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "alignments",
            "annotate",
            "XXXXXXXXXXXX",
        ],
    )

    assert result.exit_code != 0


def test_alignments_filter():

    bam = os.path.join(dir_data, "test", "Slc25A37-test_1_part0.pe.bam")
    annotations = os.path.join(
        dir_data, "test", "Slc25A37-test_1_part0.pe.annotations.tsv"
    )
    output_prefix = os.path.join(dir_test, "test", "test_filtering_cli")
    stats_prefix = os.path.join(dir_test, "stats", "test_filtering_cli")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "alignments",
            "filter",
            "capture",
            "-b",
            bam,
            "-a",
            annotations,
            "-o",
            output_prefix,
            "--stats_prefix",
            stats_prefix,
        ],
    )

    assert result.exit_code == 0

    df_slices = pd.read_csv(f"{output_prefix}.Slc25A37.slices.tsv", sep="\t")
    assert df_slices.shape[0] == 1062


def test_reporter_count():

    reporters = os.path.join(dir_data, "test", "Slc25A37_reporters.tsv.gz")
    output = os.path.join(dir_test, "test", "test_reporter_counts.tsv")
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "reporters",
            "count",
            reporters,
            "-o",
            output,
        ],
    )

    assert result.exit_code == 0


def test_reporter_storage():

    counts = os.path.join(dir_data, "test", "Slc25A37_reporter_counts.tsv.gz")
    bins = os.path.join(dir_data, "test", "genome.digest.bed.gz")
    viewpoints = os.path.join(dir_data, "test", "mm9_capture_oligos.bed")
    output_prefix = os.path.join(dir_test, "test/cli_cooler")
    output = os.path.join(dir_test, "test/cli_cooler.Slc25A37.fragments.hdf5")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "reporters",
            "store",
            "fragments",
            counts,
            "-f",
            bins,
            "-g",
            "mm9",
            "-n",
            "Slc25A37",
            "-o",
            output_prefix,
            "-c",
            viewpoints,
            "--suffix",
            "fragments",
        ],
    )

    infile = output
    output_prefix = os.path.join(dir_test, "test", "cli_cooler_binned.hdf5")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "reporters",
            "store",
            "bins",
            output,
            "-b",
            "2500",
            "-o",
            output_prefix,
            "--normalise",
            "-p",
            "4",
        ],
    )

    assert result.exit_code == 0

    runner = CliRunner()
    clrs = glob.glob(os.path.join(dir_test, "test", "cli*.hdf5"))


    result = runner.invoke(
        cli,
        [
            "reporters",
            "store",
            "merge",
            *clrs,
            '-o',
            os.path.join(dir_test, "test", "cli_cooler_merged.hdf5"),
        ],
    )

    assert result.exit_code == 0
