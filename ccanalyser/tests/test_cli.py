import os
from click.types import DateTime
import pytest
import click
from click.testing import CliRunner
import ujson
import xopen
import pysam

from ccanalyser.cli import cli

dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data")


def dicts_are_equal(d1, d2):
    return all(
        k_1 == k_2 and v_1 == v_2
        for (k_1, v_1), (k_2, v_2) in zip(d1.items(), d2.items())
    )

def get_fastq_n_records(fq):
    return sum(1 for l in pysam.FastxFile(fq))

def get_expected_duplicates():
    with open("duplicated_test_files.log") as r:
        duplicates_expected = int(
            r.readlines()[0].strip()
        )  # 4 lines per fastq

    return duplicates_expected


class DataDeduplicate:
    fq1 = os.path.join(dir_data, "duplicated_1.fastq.gz")
    fq2 = os.path.join(dir_data, "duplicated_2.fastq.gz")
    duplicates_expected = get_expected_duplicates()
    result_parsed = os.path.join(dir_test, "fq_parsed_result.json")
    result_identify = os.path.join(dir_test, "fq_duplicates_result.json")
    result_remove = os.path.join(dir_test, "fq_dedup_1.fastq")



# Test data
data_dd = DataDeduplicate


def test_cli():
    """Test checks that the cli is functional and the help option works"""
    runner = CliRunner()

    result = runner.invoke(cli, ["-h"])
    assert result.exit_code == 0

def test_genome_digest():

    test_fa = os.path.join(dir_data, 'test.fa')
    test_output = os.path.join(dir_test, 'test_digest_genome.bed')

    runner = CliRunner()
    result = runner.invoke(
            cli, ["genome-digest", test_fa, '-r', 'dpnii', '-o', test_output]
        )
    
    assert result.exit_code == 0
    
    with open(test_output) as r:
        test_bed = [l.strip().split('\t') for l in r]
    
    # Should be just the one cutsite at the moment
    assert len(test_bed) == 2

class TestFastqDeduplicate:
  

    def test_parsing(self):
        # Test parsing
        output_parsed = os.path.join(dir_test, "fq_parsed_test.json")

        runner = CliRunner()
        result = runner.invoke(
            cli, ["fastq-deduplicate", "parse", data_dd.fq1 , data_dd.fq2, "-o", output_parsed]
        )

        # Check that the script exits successfully
        assert result.exit_code == 0

        with open(output_parsed) as f:
            result_test = ujson.load(f)
        
        with open(data_dd.result_parsed) as f:
            result_correct = ujson.load(f)

        # Check that all keys and values match between the saved file and the test
        assert result_test == result_correct

    def test_identification(self):
        output_duplicates_test = os.path.join(dir_test, "fq_duplicates_test.json")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "fastq-deduplicate",
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
        
    def test_removal(self):

        # Test removal
        output_removal_prefix = os.path.join(dir_test, "fq_dedup_test")
        output_removal_test = output_removal_prefix + "_1.fastq"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "fastq-deduplicate",
                "remove",
                "-d",
                data_dd.result_identify,
                "-o",
                output_removal_prefix,
                data_dd.fq1,
                data_dd.fq2,
            ],
        )

        assert result.exit_code == 0

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

    fq1 = os.path.join(dir_data, 'digest_1.fastq.gz')
    fq2 = os.path.join(dir_data, 'digest_2.fastq.gz')
    test_output_flashed = os.path.join(dir_test, 'test_fastq_digest_flashed.fastq')
    test_output_pe = os.path.join(dir_test, 'test_fastq_digest_pe.fastq')

    runner = CliRunner()
    result = runner.invoke(
            cli, ["fastq-digest", '-m', 'pe', '-r', 'dpnii', '-o', test_output_pe, fq1, fq2]
        )
    
    assert result.exit_code == 0

    result = runner.invoke(
            cli, ["fastq-digest", '-m', 'flashed', '-r', 'dpnii', '-o', test_output_flashed, fq1]
        )
    
    assert result.exit_code == 0

