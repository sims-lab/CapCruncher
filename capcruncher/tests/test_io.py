from capcruncher.utils import PysamFakeEntry
import os
from multiprocessing import SimpleQueue
from capcruncher.tools.io import FastqReaderProcess, FastqWriterProcess, FastqWriterSplitterProcess
import itertools
import pysam

# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data")

def test_fq_reader():

    fq1 = os.path.join(dir_data, 'test', 'Slc25A37-test_1_1.fastq.gz')
    fq2 = os.path.join(dir_data, 'test', 'Slc25A37-test_1_2.fastq.gz')

    outq = SimpleQueue()
    reader = FastqReaderProcess(input_files=[fq1, fq2], outq=outq)
    reader.start()

    reads = []
    while True:
        r = outq.get()
        if not r == "END":
            reads.append(r)
        else:
            break
    
    reader.join()

    n_reads = sum(1 for read in itertools.chain.from_iterable(reads))
    assert n_reads == 1001
    
def test_fq_writer():


    fq1 = os.path.join(dir_data, 'test', 'Slc25A37-test_1_1.fastq.gz')
    fq2 = os.path.join(dir_data, 'test', 'Slc25A37-test_1_2.fastq.gz')
    fq_output = os.path.join(dir_test, 'test', 'test_fq_write.fastq')

    writeq = SimpleQueue()
    reader = FastqReaderProcess(input_files=[fq1, fq2], outq=writeq)
    writer = FastqWriterProcess(inq=writeq, output=fq_output)
    writer.start()
    reader.start()

    reader.join()
    writer.join()

    assert sum(1 for r in pysam.FastxFile(fq_output)) == 1001









