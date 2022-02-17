import logging
import pytest
import os
from click.testing import CliRunner
import glob

from capcruncher.cli import cli


@pytest.fixture(scope="module", autouse=True)
def setup_testing_dir(tmpdir_factory):
    cwd = os.getcwd()
    tmpdir = tmpdir_factory.mktemp("cli_testing")
    os.chdir(tmpdir)
    yield
    os.chdir(cwd)


@pytest.fixture(scope="session")
def testdata_dirname():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    return dirname


@pytest.fixture
def data_digestion(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "fastq_digestion")
    return data_dir


@pytest.fixture
def data_deduplication(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "fastq_deduplication")
    return data_dir


@pytest.fixture(scope="module")
def data_annotation(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "alignment_annotation")
    return data_dir


@pytest.fixture(scope="module")
def data_filter(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "alignment_filtering")
    return data_dir


@pytest.fixture(scope="module")
def data_deduplication_alignments(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "alignment_deduplication")
    return data_dir


@pytest.fixture(scope="module")
def data_reporters_count(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "reporters_count")
    return data_dir


@pytest.fixture(scope="module")
def data_reporters_store(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "reporters_store")
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
    tmp_log = os.path.join(tmpdir, "gd.log")

    result = cli_runner.invoke(
        cli, ["genome", "digest", infile, "-o", outfile, "-l", tmp_log, *flags]
    )
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


@pytest.mark.parametrize(
    "infiles,flags",
    [
        (
            ("digest_1.fastq.gz",),
            [
                "-m",
                "flashed",
                "-r",
                "dpnii",
                "--minimum_slice_length",
                "18",
            ],
        ),
        (
            ("digest_1.fastq.gz", "digest_2.fastq.gz"),
            [
                "-m",
                "pe",
                "-r",
                "dpnii",
                "--minimum_slice_length",
                "18",
            ],
        ),
        pytest.param(
            ("digest_1.fastq.gz", "digest_2.fastq.gz"),
            [
                "-m",
                "pe",
                "-r",
                "dpnii",
                "--minimum_slice_length",
                "18",
            ],
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_fastq_digest(cli_runner, data_digestion, tmpdir, infiles, flags):

    infiles = [os.path.join(data_digestion, fn) for fn in infiles]
    outfile = os.path.join(tmpdir, "digested.fastq")

    result = cli_runner.invoke(
        cli, ["fastq", "digest", *infiles, "-o", outfile, *flags]
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.mark.parametrize(
    "bam,beds,flags",
    [
        (
            "test.pe.bam",
            [
                "test_capture.bed",
            ],
            [
                "-n",
                "TEST_OVERLAP",
                "-a",
                "get",
                "-f",
                1e-9,
            ],
        ),
        (
            "test.pe.bam",
            [
                "test_capture.bed",
            ],
            [
                "-n",
                "TEST_OVERLAP",
                "-a",
                "get",
                "-f",
                1e-9,
                "--priority-chroms",
                "chr14", 
                "--prioritize-cis-slices",
            ],
        ),
    ],
)
def test_alignment_annotation(cli_runner, data_annotation, tmpdir, bam, beds, flags):

    bam = os.path.join(data_annotation, bam)
    beds = [os.path.join(data_annotation, bed) for bed in beds]
    blacklist = os.path.join(data_annotation, "test_exlcusions_corrected.bed")
    outfile = os.path.join(tmpdir, "annotated.parquet")

    result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "annotate",
            bam,
            "-b",
            *beds,
            "-o",
            outfile,
            *flags,
            "--blacklist",
            blacklist,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.mark.parametrize(
    "mode,bam,annotations,flags",
    [
        (
            "capture",
            "test.flashed.bam",
            "test.annotations.parquet",
            [],
        ),
    ],
)
def test_alignment_filter(
    cli_runner, data_filter, tmpdir, mode, bam, annotations, flags
):

    bam = os.path.join(data_filter, bam)
    annotations = os.path.join(data_filter, annotations)
    output_prefix = os.path.join(tmpdir, "filtered.parquet")
    stats_prefix = os.path.join(tmpdir, "test_filtering_cli")

    result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "filter",
            mode,
            "-b",
            bam,
            "-a",
            annotations,
            "-o",
            output_prefix,
            *flags,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(f"{output_prefix}.slices.parquet")


@pytest.mark.parametrize(
    "fragments,output",
    [
        (
            "fragments.flashed.parquet",
            "dup.flashed.pkl",
        ),
        (
            "fragments.pe.parquet",
            "dup.pe.pkl",
        ),
    ],
)
def test_alignment_deduplicate_fragments(
    cli_runner, data_deduplication_alignments, tmpdir, fragments, output
):

    fragments = os.path.join(data_deduplication_alignments, fragments)
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "deduplicate",
            "identify",
            fragments,
            "-o",
            output,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(output)


@pytest.mark.parametrize(
    "slices,duplicates,output",
    [
        (
            "slices.flashed.parquet",
            "duplicates.flashed.pkl",
            "deduplicated.flashed.parquet",
        ),
        (
            "slices.pe.parquet",
            "duplicates.pe.pkl",
            "deduplicated.pe.parquet",
        ),
    ],
)
def test_alignment_deduplicate_slices(
    cli_runner, data_deduplication_alignments, tmpdir, slices, duplicates, output
):

    slices = os.path.join(data_deduplication_alignments, slices)
    duplicates = os.path.join(data_deduplication_alignments, duplicates)
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "alignments",
            "deduplicate",
            "remove",
            slices,
            "-d",
            duplicates,
            "-o",
            output,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(output)


@pytest.mark.parametrize(
    "slices,bins,viewpoints,output,flags",
    [
        (
            "slices.flashed.parquet",
            "bins.bed.gz",
            "viewpoints.bed",
            "counts.hdf5",
            ["--cooler-output"],
        ),
    ],
)
def test_reporters_count(
    cli_runner,
    data_deduplication_alignments,
    data_reporters_count,
    tmpdir,
    slices,
    bins,
    viewpoints,
    output,
    flags,
):

    slices = os.path.join(data_deduplication_alignments, slices)
    bins = os.path.join(data_reporters_count, bins)
    viewpoints = os.path.join(data_reporters_count, viewpoints)
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "reporters",
            "count",
            slices,
            "-o",
            output,
            "-f",
            bins,
            "-v",
            viewpoints,
            *flags,
        ],
    )

    logging.debug(result.stdout)
    logging.debug(result.exception)
    assert result.exit_code == 0
    assert os.path.exists(output)


@pytest.mark.parametrize(
    "cooler_fn,bin_size,output,flags",
    [
        ("SAMPLE-A_REP1.hdf5", int(1e5), "binned.hdf5", []),
    ],
)
def test_reporters_store_binned(
    cli_runner,
    data_reporters_store,
    tmpdir,
    cooler_fn,
    bin_size,
    output,
    flags,
):

    clr = os.path.join(data_reporters_store, cooler_fn)
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "reporters",
            "store",
            "bins",
            clr,
            "-o",
            output,
            "-b",
            bin_size,
            *flags,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(output)


@pytest.mark.parametrize(
    "cooler_fn,output_prefix,outfile,flags",
    [
        ("SAMPLE-A_REP1.hdf5", "test", "test.Slc25A37.bedgraph", []),
        ("SAMPLE-A_REP1.hdf5", "test", "test.Slc25A37.bigWig", ["-f", "bigwig"]),
    ],
)
def test_reporters_pileup(
    cli_runner,
    data_reporters_store,
    tmpdir,
    cooler_fn,
    output_prefix,
    outfile,
    flags,
):

    clr = os.path.join(data_reporters_store, cooler_fn)
    output = os.path.join(tmpdir, output_prefix)
    outfile = os.path.join(tmpdir, outfile)

    result = cli_runner.invoke(
        cli,
        [
            "reporters",
            "pileup",
            clr,
            "-o",
            output,
            *flags,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)
