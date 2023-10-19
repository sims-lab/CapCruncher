import pandas as pd
import pysam
import pytest
import os


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "fastq_digestion")
    return data_dir


def count_fragments(fq):
    fragments = set()
    with pysam.FastxFile(fq) as fq_file:
        for record in fq_file:
            parent_id = record.name.split("|")[0]
            fragments.add(parent_id)
    return len(fragments)


@pytest.mark.parametrize(
    "fastq_files,enzyme,mode,n_reads_raw,n_reads_filt",
    [
        # (("test.fq.gz.extendedFrags.fastq.gz",), "dpnii", "flashed", 671471, 478386),
        (("digest_1.fastq.gz",), "dpnii", "flashed", 1512, 876),
        (("digest_1.fastq.gz", "digest_2.fastq.gz"), "dpnii", "pe", 1512, 1512),
        pytest.param(
            ("digest_1.fastq.gz", "digest_2.fastq.gz"),
            "dpnii",
            "flashed",
            1512,
            876,
            marks=pytest.mark.xfail(strict=True),
            id="flashed_should_fail",
        ),
        pytest.param(
            ("digest_1.fastq.gz",),
            "dpnii",
            "pe",
            1512,
            876,
            marks=pytest.mark.xfail(strict=True),
            id="pe_should_fail",
        ),
    ],
)
def test_digest_fastq(
    data_path, tmpdir, fastq_files, enzyme, mode, n_reads_raw, n_reads_filt
):
    from capcruncher.cli.fastq_digest import digest

    infiles = [os.path.join(data_path, fn) for fn in fastq_files]
    outfile = os.path.join(tmpdir, "out.fq")
    stats_prefix = os.path.join(tmpdir, "stats")

    stats = digest(
        infiles,
        enzyme,
        mode=mode,
        output_file=outfile,
        stats_prefix=stats_prefix,
    )

    assert (
        stats["stats_read_level"].to_pandas()["number_of_read_pairs_unfiltered"].iloc[0]
        == n_reads_raw
    )
    assert (
        stats["stats_read_level"].to_pandas()["number_of_read_pairs_filtered"].iloc[0]
        == n_reads_filt
    )
    assert count_fragments(outfile) == n_reads_filt


@pytest.fixture(scope="module")
def fasta():
    import pathlib

    fa = (
        pathlib.Path(__file__).parent / "data" / "data_for_pipeline_run" / "chr14.fa.gz"
    )
    return str(fa)


@pytest.mark.parametrize(
    "enzyme,n_records_expected",
    [
        ("dpnii", 2),
    ],
)
def test_digest_genome(fasta, tmpdir, enzyme, n_records_expected):
    from capcruncher.cli.genome_digest import digest

    infile = fasta
    outfile = os.path.join(tmpdir, "digested.bed")

    digest(input_fasta=infile, recognition_site=enzyme, output_file=outfile)

    assert os.path.exists(outfile)
