import multiprocessing
import pandas as pd
import pysam
import pytest
import os

from click.testing import CliRunner
from capcruncher.tools.digest import DigestedRead, DigestedChrom, ReadDigestionProcess
from capcruncher.utils import MockFastaRecord, MockFastqRecord
from multiprocessing import SimpleQueue


def make_fake_read(sequence):
    return MockFastqRecord("read", sequence, "".join(["A" for _ in sequence]))


def make_fake_fasta(sequence):
    return MockFastaRecord("chrom", sequence)


@pytest.mark.parametrize(
    "sequence,allow_undigested,cutsite,min_slice_len,slice_lengths,n_valid",
    [
        (
            "".join(["A" * 20, "G" * 20, "GATC", "T" * 20, "C" * 20]),
            False,
            "GATC",
            18,
            [0, 40, 84],
            2,
        ),
        (
            "".join(["A" * 20, "G" * 20, "T" * 20, "C" * 20]),
            True,
            "GATC",
            18,
            [0, 80],
            1,
        ),
        (
            "".join(["A" * 20, "G" * 20, "T" * 20, "C" * 20]),
            False,
            "GATC",
            18,
            [0, 80],
            0,
        ),
        (
            "".join(["A" * 5, "GATC", "G" * 20, "T" * 20, "C" * 20]),
            False,
            "GATC",
            18,
            [0, 5, 69],
            1,
        ),
        (
            "".join(["A" * 5, "GATC", "G" * 20, "T" * 20, "C" * 20]),
            True,
            "GATC",
            18,
            [0, 5, 69],
            1,
        ),
    ],
)
def test_digest_reads(
    sequence, allow_undigested, cutsite, min_slice_len, slice_lengths, n_valid
):

    read = make_fake_read(sequence)

    dr = DigestedRead(
        read=read,
        cutsite=cutsite,
        allow_undigested=allow_undigested,
        min_slice_length=min_slice_len,
    )
    dr_indexes = dr.get_recognition_site_indexes()

    assert dr_indexes == slice_lengths
    assert dr.slices_filtered == n_valid


@pytest.mark.parametrize(
    "sequence,allow_undigested,cutsite,min_slice_len,slice_lengths,n_valid",
    [
        (
            "".join(["A" * 20, "G" * 20, "GATC", "T" * 20, "C" * 20]),
            False,
            "GATC",
            18,
            [40, 40],
            2,
        ),
        (
            "".join(["A" * 20, "G" * 20, "T" * 20, "C" * 20]),
            True,
            "GATC",
            18,
            [80],
            1,
        ),
        (
            "".join(["A" * 20, "G" * 20, "T" * 20, "C" * 20]),
            False,
            "GATC",
            18,
            [],
            0,
        ),
        (
            "".join(["A" * 5, "GATC", "G" * 20, "T" * 20, "C" * 20]),
            False,
            "GATC",
            18,
            [60],
            1,
        ),
        (
            "".join(["A" * 5, "GATC", "G" * 20, "T" * 20, "C" * 20]),
            True,
            "GATC",
            18,
            [60],
            1,
        ),
    ],
)
def test_digest_reads_process(
    sequence, allow_undigested, cutsite, min_slice_len, slice_lengths, n_valid
):

    read = make_fake_read(sequence)

    inq = multiprocessing.Queue()
    outq = multiprocessing.Queue()

    dp = ReadDigestionProcess(
        inq,
        outq,
        cutsite=cutsite,
        allow_undigested=allow_undigested,
        min_slice_length=min_slice_len,
    )

    dp.start()
    inq.put(
        [
            (read,),
        ]
    )
    inq.put(None)
    digested = outq.get()
    dp.join()

    digested_lines_sequences = digested.split("\n")[1::4]
    assert [len(s) for s in digested_lines_sequences] == slice_lengths
    assert len(digested_lines_sequences) == n_valid


@pytest.mark.parametrize(
    "sequence,cutsite,indexes,coordinates",
    [
        (
            "".join(["A" * 20, "G" * 20, "GATC", "T" * 20, "C" * 20]),
            "GATC",
            [0, 40, 84],
            "chrom\t0\t40\t0\nchrom\t44\t84\t1\n",
        ),
    ],
)
def test_digest_chrom(sequence, cutsite, indexes, coordinates):

    record = make_fake_fasta(sequence)
    dr = DigestedChrom(chrom=record, cutsite=cutsite)
    assert dr.fragment_indexes == indexes
    assert "".join(dr.fragments) == coordinates


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
        #(("test.fq.gz.extendedFrags.fastq.gz",), "dpnii", "flashed", 671471, 478386),
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
        n_cores=3,
    )

    test_n_reads_raw = stats.query("(stat_type == 'unfiltered') and (read_number < 2)")["stat"].values[0]
    test_n_reads_filt = stats.query("(stat_type == 'filtered') and (read_number < 2)")["stat"].values[0]
    hist_raw = pd.read_csv(f"{stats_prefix}.digestion.unfiltered.histogram.csv")
    hist_filt = pd.read_csv(f"{stats_prefix}.digestion.filtered.histogram.csv")

    assert test_n_reads_raw == n_reads_raw
    assert test_n_reads_filt == n_reads_filt
    assert count_fragments(outfile) == n_reads_filt
    assert hist_raw.query("read_number < 2")["count"].sum() == n_reads_raw 
    assert hist_filt.query("read_number < 2")["count"].sum() == n_reads_filt   




@pytest.mark.parametrize(
    "fasta,enzyme,n_records_expected",
    [
        ("chrom_to_digest.fa", "dpnii", 2),
    ],
)
def test_digest_genome(
    data_path, tmpdir, fasta, enzyme, n_records_expected
):

    from capcruncher.cli.genome_digest import digest

    infile = os.path.join(data_path, fasta)
    outfile = os.path.join(tmpdir, "digested.bed")

    digest(input_fasta=infile, recognition_site=enzyme, output_file=outfile)

    assert os.path.exists(outfile)


