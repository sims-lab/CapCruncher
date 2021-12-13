import multiprocessing
import queue
from capcruncher.utils import PysamFakeEntry
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
multiprocessing.set_start_method("fork")
import gzip

# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data", "test", "io")


@pytest.mark.parametrize(
    "fastq_files,n_records",
    [
        (
            (
                os.path.join(dir_data, "Slc25A37-test_1_1.fastq.gz"),
                os.path.join(dir_data, "Slc25A37-test_1_2.fastq.gz"),
            ),
            1001,
        ),
        (
            (os.path.join(dir_data, "Slc25A37-test_1_1.fastq.gz"),),
            1001,
        ),
    ],
)
def test_fq_reader(fastq_files, n_records):

    outq = multiprocessing.Queue()
    reader = FastqReaderProcess(input_files=fastq_files, outq=outq)
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
                os.path.join(dir_data, "Slc25A37-test_1_1.fastq.gz"),
                os.path.join(dir_data, "Slc25A37-test_1_2.fastq.gz"),
            ),
            (
                os.path.join(dir_test, "written_Slc25A37-test_1_1.fastq.gz"),
                os.path.join(dir_test, "written_Slc25A37-test_1_2.fastq.gz"),
            ),
            1001,
        ),
        (
            (os.path.join(dir_data, "Slc25A37-test_1_1.fastq.gz"),),
            (os.path.join(dir_test, "written_Slc25A37-test_single.fastq.gz"),),
            1001,
        ),
    ],
)
def test_fq_writer(in_files, out_files, n_records_expected):

    outq = multiprocessing.Queue()
    reader = FastqReaderProcess(input_files=in_files, outq=outq)
    writer = FastqWriterProcess(inq=outq, output=out_files)

    reader.start()
    writer.start()

    reader.join()
    writer.join()
 
    n_records_test = len([rec.name for rec in pysam.FastxFile(out_files[0])])  
    assert n_records_test == n_records_expected
    

