import multiprocessing
from posixpath import dirname
import queue
from capcruncher.utils import MockFastqRecord
import os
from multiprocessing import SimpleQueue
from capcruncher.tools.io import (
    FastqReaderProcess,
    FastqWriterProcess,
    FastqWriterSplitterProcess,
)
import itertools
import pysam
import pytest
import xopen
import gzip

# Pre-run setup


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "io")
    return data_dir


@pytest.mark.parametrize(
    "fastq_files,n_records",
    [
        (
            (
                "Slc25A37-test_1_1.fastq.gz",
                "Slc25A37-test_1_2.fastq.gz",
            ),
            1001,
        ),
        (
            ("Slc25A37-test_1_1.fastq.gz",),
            1001,
        ),
    ],
)
def test_fq_reader(data_path, fastq_files, n_records):

    outq = multiprocessing.Queue()
    reader = FastqReaderProcess(
        input_files=[os.path.join(data_path, fn) for fn in fastq_files], outq=outq
    )
    reader.start()

    reads = []
    while True:

        try:
            r = outq.get(block=True, timeout=0.01)
        except queue.Empty:
            continue

        if r:
            reads.append(r)
        else:
            break

    reader.join()

    n_reads = len(list(itertools.chain.from_iterable(reads)))
    assert n_reads == n_records


@pytest.mark.parametrize(
    "in_files,out_files,n_records_expected",
    [
        (
            (
                "Slc25A37-test_1_1.fastq.gz",
                "Slc25A37-test_1_2.fastq.gz",
            ),
            (
                "written_Slc25A37-test_1_1.fastq.gz",
                "written_Slc25A37-test_1_2.fastq.gz",
            ),
            1001,
        ),
        (
            ("Slc25A37-test_1_1.fastq.gz",),
            ("written_Slc25A37-test_single.fastq.gz",),
            1001,
        ),
    ],
)
def test_fq_writer(data_path, tmpdir, in_files, out_files, n_records_expected):

    outq = multiprocessing.Queue()

    in_files = [os.path.join(data_path, fn) for fn in in_files]
    out_files = [os.path.join(tmpdir, fn) for fn in out_files]

    reader = FastqReaderProcess(input_files=in_files, outq=outq)
    writer = FastqWriterProcess(inq=outq, output=out_files)

    reader.start()
    writer.start()

    reader.join()
    writer.join()

    n_records_test = len([rec.name for rec in pysam.FastxFile(out_files[0])])
    assert n_records_test == n_records_expected
