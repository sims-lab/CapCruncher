from loguru import logger
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
def data_pipeline(testdata_dirname):
    data_dir = os.path.join(testdata_dirname, "data", "data_for_pipeline_run")
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
            "chr14.fa.gz",
            ["-r", "dpnii", "--sort"],
        ),
        pytest.param(
            "chr14.fa.gz",
            [
                "-r",
                "dpn",
            ],
            marks=pytest.mark.xfail,
        ),
    ],
)
def test_genome_digest(cli_runner, data_pipeline, tmpdir, infile, flags):
    infile = os.path.join(data_pipeline, infile)
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
            "out",
            [],
        )
    ],
)
def test_deduplicate_fastq(
    cli_runner, data_deduplication, tmpdir, infiles, outfile, flags
):
    infiles = [os.path.join(data_deduplication, fn) for fn in infiles]
    outfile = os.path.join(tmpdir, outfile)

    result = cli_runner.invoke(
        cli,
        [
            "fastq",
            "deduplicate",
            "-1",
            infiles[0],
            "-2",
            infiles[1],
            "-o",
            outfile,
            *flags,
        ],
    )

    assert result.exit_code == 0

    fastq_deduped = glob.glob(f"{tmpdir}/*.fastq*")
    assert len(fastq_deduped) == len(infiles)


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
        (
            "tri",
            "test.flashed.bam",
            "test.annotations.parquet",
            [],
        ),
        (
            "tiled",
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
            "--sample-name",
            "test",
            *flags,
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0
    assert os.path.exists(f"{output_prefix}.slices.parquet")


@pytest.mark.parametrize(
    "slices,read_type,output",
    [
        (
            "slices.flashed.parquet",
            "flashed",
            "deduplicated.flashed.parquet",
        ),
        (
            "slices.pe.parquet",
            "pe",
            "deduplicated.pe.parquet",
        ),
    ],
)
def test_interactions_deduplicate(
    cli_runner, data_deduplication_alignments, tmpdir, slices, read_type, output
):
    slices = os.path.join(data_deduplication_alignments, slices)
    output = os.path.join(tmpdir, output)

    """
    slices: os.PathLike,
    output: os.PathLike,
    read_type: str = "flashed",
    sample_name: str = "sampleX",
    statistics: os.PathLike = "deduplication_stats.json",
    """

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
            "deduplicate",
            slices,
            "--read-type",
            read_type,
            "-o",
            output,
            "--sample-name",
            "test",
            "--statistics",
            "test_deduplication_cli.json",
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(output)


@pytest.mark.parametrize(
    "slices,bins,viewpoints,output,flags",
    [
        ("slices.flashed.parquet", "bins.bed.gz", "viewpoints.bed", "counts.hdf5", []),
    ],
)
def test_reporters_count(
    cli_runner,
    data_deduplication_alignments,
    data_reporters_count,
    data_pipeline,
    tmpdir,
    slices,
    bins,
    viewpoints,
    output,
    flags,
):
    slices = os.path.join(data_deduplication_alignments, slices)
    bins = os.path.join(data_reporters_count, bins)
    viewpoints = os.path.join(data_pipeline, "mm9_capture_viewpoints_Slc25A37.bed")
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
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

    logger.debug(result.stdout)
    logger.debug(result.exception)
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
            "interactions",
            "fragments-to-bins",
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
    "infiles,viewpoint,output,flags",
    [
        (["SAMPLE-A_REP1.hdf5", "SAMPLE-A_REP1.hdf5"], "Slc25A37", "merged.hdf5", []),
    ],
)
def test_reporters_store_merge(
    cli_runner,
    data_reporters_store,
    tmpdir,
    infiles,
    viewpoint,
    output,
    flags,
):
    import cooler

    hdf5_files = [os.path.join(data_reporters_store, fn) for fn in infiles]
    output = os.path.join(tmpdir, output)

    result = cli_runner.invoke(
        cli,
        [
            "interactions",
            "merge",
            *hdf5_files,
            "-o",
            output,
            *flags,
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(output)
    assert (
        cooler.Cooler(f"{output}::{viewpoint}").pixels()[:]["count"].sum()
        == cooler.Cooler(f"{hdf5_files[0]}::{viewpoint}").pixels()[:]["count"].sum() * 2
    )


@pytest.mark.parametrize(
    "cooler_fn,output_prefix,outfile,flags",
    [
        ("SAMPLE-A_REP1.hdf5", "test", "test_Slc25A37.bedgraph", []),
        ("SAMPLE-A_REP1.hdf5", "test", "test_Slc25A37.bigWig", ["-f", "bigwig"]),
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
            "interactions",
            "pileup",
            clr,
            "-o",
            output,
            *flags,
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(outfile)


@pytest.fixture
def fragments_file(tmpdir):
    fragments = os.path.join(tmpdir, "fragments.bed")
    with open(fragments, "w") as file:
        file.write("chr1\t100\t200\tfragment1\n")
        file.write("chr2\t300\t400\tfragment2\n")
    return fragments


@pytest.fixture
def viewpoints_file(tmpdir):
    viewpoints = os.path.join(tmpdir, "viewpoints.bed")
    with open(viewpoints, "w") as file:
        file.write("chr1\t150\t160\tviewpoint1\n")
        file.write("chr2\t350\t360\tviewpoint2\n")
    return viewpoints


def test_make_chicago_maps(cli_runner, tmpdir, fragments_file, viewpoints_file):
    outputdir = str(tmpdir)

    result = cli_runner.invoke(
        cli,
        [
            "utilities",
            "make-chicago-maps",
            "--fragments",
            fragments_file,
            "--viewpoints",
            viewpoints_file,
            "-o",
            outputdir,
        ],
    )

    assert result.exit_code == 0

    # Check if the renamed fragments file exists
    fragments_new = os.path.join(outputdir, "fragments.rmap")
    assert os.path.exists(fragments_new)

    # Check if the baitmap file exists and has the correct content
    baitmap_file = os.path.join(outputdir, "viewpoints.baitmap")
    assert os.path.exists(baitmap_file)

    with open(baitmap_file, "r") as file:
        content = file.read()
        assert "chr1\t100\t200\tfragment1\tviewpoint" in content
