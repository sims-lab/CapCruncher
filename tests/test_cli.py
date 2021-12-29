import pytest
import os
from click.testing import CliRunner
import glob

from capcruncher.cli import cli


@pytest.fixture
def data_digestion():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "fastq_digestion")
    return data_dir


@pytest.fixture
def data_deduplication():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "fastq_deduplication")
    return data_dir


@pytest.fixture(scope="module")
def cli_runner():
    return CliRunner()


def test_cli_runs(cli_runner):
    """Test checks that the cli is functional and the help option works"""

    result = cli_runner.invoke(cli, ["--help"])
    assert result.exit_code == 0


@pytest.mark.parametrize(
    "infile,flags",
    [
        (
            "chrom_to_digest.fa",
            ["-r", "dpnii", "--sort"],
        ),
        pytest.param(
            "chrom_to_digest.fa",
            [
                "-r",
                "dpn",
            ],
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_genome_digest(cli_runner, data_digestion, tmpdir, infile, flags):

    infile = os.path.join(data_digestion, infile)
    outfile = os.path.join(tmpdir, "digested.bed")

    result = cli_runner.invoke(cli, ["genome", "digest", infile, "-o", outfile, *flags])
    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.mark.parametrize(
    "infiles,outfile,flags",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "parsed.json",
            ["--read_buffer", "100000"],
        ),
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "parsed.pkl",
            ["--read_buffer", "100000"],
        ),
    ],
)
def test_deduplicate_parse(
    cli_runner, data_deduplication, tmpdir, infiles, outfile, flags
):

    infiles = [os.path.join(data_deduplication, fn) for fn in infiles]
    outfile = os.path.join(tmpdir, outfile)
    result = cli_runner.invoke(
        cli, ["fastq", "deduplicate", "parse", *infiles, "-o", outfile, *flags]
    )

    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.mark.parametrize(
    "infiles,outfile,flags",
    [
        (("parsed.json",), "deduplicated.json", []),
        (("parsed.pickle",), "deduplicated.pickle", []),
        (("parsed.json", "parsed.pickle"), "deduplicated.pickle", []),
    ],
)
def test_deduplicate_identify(
    cli_runner, data_deduplication, tmpdir, infiles, outfile, flags
):

    infiles = [os.path.join(data_deduplication, fn) for fn in infiles]
    outfile = os.path.join(tmpdir, outfile)
    result = cli_runner.invoke(
        cli, ["fastq", "deduplicate", "identify", *infiles, "-o", outfile, *flags]
    )

    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.mark.parametrize(
    "infiles,ids,outfile,flags",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "identified.pkl",
            "out",
            [],
        )
    ],
)
def test_deduplicate_remove(
    cli_runner, data_deduplication, tmpdir, infiles, ids, outfile, flags
):

    infiles = [os.path.join(data_deduplication, fn) for fn in infiles]
    ids = os.path.join(data_deduplication, ids)
    outfile = os.path.join(tmpdir, outfile)
    result = cli_runner.invoke(
        cli,
        ["fastq", "deduplicate", "remove", *infiles, "-d", ids, "-o", outfile, *flags],
    )

    assert result.exit_code == 0
    assert len(glob.glob(f"{tmpdir}/*.fastq*")) == 2
