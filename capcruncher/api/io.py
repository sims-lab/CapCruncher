import glob
from loguru import logger
import multiprocessing
import os
import pathlib
import queue
import random
import string
import traceback
from typing import Literal, Union

import pandas as pd
import pysam
import tqdm
from capcruncher.utils import get_timing
from pysam import FastxFile
from xopen import xopen
import xxhash
from collections import defaultdict, namedtuple


class FastqReaderProcess(multiprocessing.Process):
    """Reads fastq file(s) in chunks and places them on a queue.

    Attributes:
     input_file: Input fastq files.
     outq: Output queue for chunked reads/read pairs.
     statq: (Not currently used) Queue for read statistics if required.
     read_buffer: Number of reads to process before placing them on outq
     read_counter: (Not currently used) Can be used to sync between multiple readers.
     n_subproceses: Number of processes running concurrently. Used to make sure enough termination signals are used.

    """

    def __init__(
        self,
        input_files: Union[str, list],
        outq: multiprocessing.Queue,
        read_buffer: int = 100000,
    ) -> None:
        # Input variables
        self.input_files = input_files
        self._multifile = self._is_multifile(input_files)

        if self._multifile:
            self._input_files_pysam = [FastxFile(f) for f in self.input_files]
        else:
            self._input_files_pysam = [
                FastxFile(self.input_files),
            ]

        # Multiprocessing variables
        self.outq = outq

        # Reader variables
        self.read_buffer = read_buffer

        super(FastqReaderProcess, self).__init__()

    def _is_multifile(self, files):
        if not isinstance(files, (str, pathlib.Path)):
            return True
        elif isinstance(files, (list, tuple)) and len(files > 1):
            return True
        else:
            return False

    def run(self):
        """Performs reading and chunking of fastq file(s)."""

        try:
            buffer = []
            rc = 0
            for read_counter, read in enumerate(zip(*self._input_files_pysam)):
                # print(f"read_counter: {read_counter}, read: {read}, read_buffer: {self.read_buffer}")
                buffer.append(read)
                if read_counter % self.read_buffer == 0 and not read_counter == 0:
                    self.outq.put(buffer.copy())
                    buffer.clear()
                    logger.info(f"{read_counter} reads parsed (batch)")
                    rc = read_counter
                else:
                    rc = read_counter

            self.outq.put(buffer)  # Deal with remainder
            self.outq.put_nowait(None)  # Poison pill to terminate queue
            logger.info(f"{rc} reads parsed (final)")

        except Exception as e:
            logger.info(f"Reader failed with exception: {e}")
            raise

        finally:
            for fh in self._input_files_pysam:
                fh.close()


class FastqReadFormatterProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.SimpleQueue,
        outq: multiprocessing.SimpleQueue,
        formatting: list = None,
    ) -> None:
        self.inq = inq
        self.outq = outq
        self.formatting = (
            [
                self._format_as_str,
            ]
            if not formatting
            else formatting
        )

        super(FastqReadFormatterProcess, self).__init__()

    def _format_as_str(self, reads):
        # [(r1, r2), (r1, r2)] -> [r1 combined string, r2 combined string]
        return ["\n".join([str(rn) for rn in r]) for r in zip(*reads)]

    def run(self):
        try:
            reads = self.inq.get()

            while not reads == "END":
                for formatting_to_apply in self.formatting:
                    reads = formatting_to_apply(reads)

                self.outq.put(reads)
                reads = self.inq.get()

            self.outq.put("END")

        except Exception:
            traceback.format_exc()
            self.outq.put("END")


class FastqWriterSplitterProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        output_prefix: Union[str, list],
        paired_output: bool = False,
        gzip=False,
        compression_level: int = 3,
        compression_threads: int = 8,
        n_subprocesses: int = 1,
        n_workers_terminated: int = 0,
        n_files_written: int = 0,
    ):
        self.inq = inq
        self.output_prefix = output_prefix
        self.paired_output = paired_output

        self.gzip = gzip
        self.compression_level = compression_level
        self.compression_threads = compression_threads

        self.n_subprocesses = n_subprocesses
        self.n_workers_terminated = n_workers_terminated
        self.n_files_written = n_files_written

        super(FastqWriterSplitterProcess, self).__init__()

    def _get_file_handles(self):
        if not self.paired_output:
            fnames = [
                f'{self.output_prefix}_part{self.n_files_written}.fastq{".gz" if self.gzip else ""}',
            ]
        else:
            fnames = [
                f'{self.output_prefix}_part{self.n_files_written}_{i+1}.fastq{".gz" if self.gzip else ""}'
                for i in range(2)
            ]

        return [
            xopen(
                fn,
                "w",
                compresslevel=self.compression_level,
                threads=self.compression_threads,
            )
            for fn in fnames
        ]

    def run(self):
        try:
            reads = self.inq.get()
            is_string_input = True if isinstance(reads[0], str) else False

            while self.n_workers_terminated < self.n_subprocesses:
                if reads == "END":
                    self.n_workers_terminated += 1
                    continue

                elif is_string_input:
                    for fh, read in zip(self._get_file_handles(), reads):
                        fh.write(read)
                        fh.close()

                else:
                    reads_str = [
                        "\n".join([str(r) for r in read_glob])
                        for read_glob in zip(*reads)
                    ]

                    for fh, read_set in zip(self._get_file_handles(), reads_str):
                        fh.write((read_set + "\n"))
                        fh.close()

                reads = self.inq.get()
                self.n_files_written += 1

        except Exception:
            traceback.format_exc()


CCAlignment = namedtuple(
    "CCAlignment",
    field_names=[
        "slice_id",
        "slice_name",
        "parent_id",
        "parent_read",
        "pe",
        "slice",
        "uid",
        "mapped",
        "multimapped",
        "chrom",
        "start",
        "end",
        "coordinates",
    ],
)


def parse_alignment(aln: pysam.AlignmentFile) -> CCAlignment:
    """Parses reads from a bam file into a list.

    Extracts:
     -read name
     -parent reads
     -flashed status
     -slice number
     -mapped status
     -multimapping status
     -chromosome number (e.g. chr10)
     -start (e.g. 1000)
     -end (e.g. 2000)
     -coords e.g. (chr10:1000-2000)


    Args:
     aln: pysam.AlignmentFile.
    Returns:
     list: Containing the attributes extracted.

    """

    import numpy as np

    slice_name = aln.query_name
    parent_read, pe, slice_number, uid = slice_name.split("|")
    parent_id = xxhash.xxh3_64_intdigest(parent_read, seed=42)
    slice_id = xxhash.xxh3_64_intdigest(slice_name, seed=42)
    ref_name = aln.reference_name
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    # Check if read mapped
    if aln.is_unmapped:
        mapped = 0
        multimapped = 0
        ref_name = ""
        ref_start = 0
        ref_end = 0
        coords = ""
    else:
        mapped = 1
        coords = f"{ref_name}:{ref_start}-{ref_end}"
        # Check if multimapped
        if aln.is_secondary:
            multimapped = 1
        else:
            multimapped = 0

    return CCAlignment(
        slice_id=slice_id,
        slice_name=slice_name,
        parent_id=parent_id,
        parent_read=parent_read,
        pe=pe.lower(),
        slice=int(slice_number),
        uid=int(uid),
        mapped=mapped,
        multimapped=multimapped,
        chrom=ref_name,
        start=int(ref_start),
        end=int(ref_end),
        coordinates=coords,
    )


def parse_bam(bam: Union[str, pathlib.Path]) -> pd.DataFrame:
    """Uses parse_alignment function convert bam file to a dataframe.

    Extracts:
     -'slice_name'
     -'parent_read'
     -'pe'
     -'slice'
     -'mapped'
     -'multimapped'
     -'chrom'
     -'start'
     -'end'
     -'coordinates'

    Args:
        bam: Path to bam file.

    Returns:
     pd.Dataframe: DataFrame with the columns listed above.

    """

    import numpy as np

    # Load reads into dataframe
    logger.info("Parsing BAM file")
    df_bam = pd.DataFrame(
        [
            parse_alignment(aln)
            for aln in pysam.AlignmentFile(bam, "rb").fetch(until_eof=True)
        ],
    )
    df_bam["bam"] = os.path.basename(bam)

    # Perform dtype conversions
    logger.info("Converting dtypes")
    df_bam["chrom"] = df_bam["chrom"].astype("category")
    pe_category = pd.CategoricalDtype(["flashed", "pe"])
    df_bam["pe"] = df_bam["pe"].astype(
        pe_category
    )  # Only the one type present so need to include both
    df_bam["coordinates"] = df_bam["coordinates"].astype("category")
    df_bam["parent_read"] = df_bam["parent_read"].astype("category")
    df_bam["slice"] = df_bam["slice"].astype(np.int8)
    df_bam["uid"] = df_bam["uid"].astype(np.int8)
    df_bam["multimapped"] = df_bam["multimapped"].astype(bool)
    df_bam["mapped"] = df_bam["mapped"].astype(bool)
    df_bam["bam"] = df_bam["bam"].astype("category")

    logger.info("Finished parsing BAM file")
    return df_bam


def bam_to_parquet(
    bam: Union[str, pathlib.Path], output: Union[str, pathlib.Path]
) -> Union[str, pathlib.Path]:
    """Converts bam file to parquet file.

    Args:
     bam: Path to bam file.
     output: Path to output parquet file.

    """
    df_bam = parse_bam(bam)
    df_bam.to_parquet(output)

    return output
