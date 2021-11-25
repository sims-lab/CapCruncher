import pathlib
import multiprocessing
from typing import Dict, Union
import traceback

import pandas as pd
from pysam import FastxFile
from xopen import xopen
from capcruncher.utils import get_timing, hash_column, get_fragment_combinations
import pysam
import numpy as np
import queue
import os


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
        read_counter: multiprocessing.Manager().Value = None,
        n_subprocesses: int = 1,
        statq: multiprocessing.Queue = None,
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
        self.statq = statq
        self.n_subprocesses = n_subprocesses

        # Reader variables
        self.read_buffer = read_buffer
        self.read_counter = read_counter

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
            for read_counter, read in enumerate(zip(*self._input_files_pysam)):

                buffer.append(read)

                if read_counter % self.read_buffer == 0 and not read_counter == 0:
                    self.outq.put(buffer)
                    buffer = []
                    print(f"Read {read_counter} reads")

            self.outq.put(buffer)
            print(f"Number of reads processed: {read_counter + 1}")

            for i in range(self.n_subprocesses):
                self.outq.put("END")

            # Deal with number of reads that have been read
            if self.read_counter:
                self.read_counter.value = read_counter
            elif self.statq:
                self.statq.put({"reads_total": read_counter})

        except Exception as e:
            traceback.format_exc()
            self.outq.put("END")


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

        except Exception as e:
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

        except Exception as e:
            traceback.format_exc()


class FastqWriterProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        output: Union[str, list],
        compression_level: int = 5,
        n_subprocesses: int = 1,
    ):

        super(FastqWriterProcess, self).__init__()

        self.inq = inq
        self.output = output
        self.compression_level = compression_level
        self.n_workers_terminated = 0
        self.n_subprocesses = n_subprocesses
        self.file_handles = self._get_filehandles()
        self.name = "FastqWriter"

    def _get_filehandles(self):
        if isinstance(self.output, str):
            return [
                xopen(
                    self.output, "w", compresslevel=self.compression_level, threads=2
                ),
            ]
        elif isinstance(self.output, (list, tuple, pd.Series)):
            return [
                xopen(fn, "w", compresslevel=self.compression_level, threads=2)
                for fn in self.output
            ]

    def _inputs_match_number_of_handles(self, reads):
        if len(reads[0]) == len(self.output):
            return True

    def run(self):

        try:

            reads = self.inq.get()
            is_string_input = True if isinstance(reads, str) else False

            counter = 0
            while self.n_workers_terminated < self.n_subprocesses:

                if reads == "END":
                    self.n_workers_terminated += 1
                    continue

                elif is_string_input:
                    for fh in self.file_handles:
                        fh.write(reads)

                else:
                    reads_str = [
                        "\n".join([str(r) for r in read_glob])
                        for read_glob in zip(*reads)
                    ]

                    for fh, read_set in zip(self.file_handles, reads_str):
                        fh.write((read_set + "\n"))

                reads = self.inq.get()

            for fh in self.file_handles:
                fh.close()

        except Exception as e:
            traceback.format_exc()
            self.outq.put("END")


def parse_alignment(aln):
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

    slice_name = aln.query_name
    parent_read, pe, slice_number, uid = slice_name.split("|")
    ref_name = aln.reference_name
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    # Check if read mapped
    if aln.is_unmapped:
        mapped = 0
        multimapped = 0
        ref_name = "unmapped"
        ref_start = ""
        ref_end = ""
        coords = ""
    else:
        mapped = 1
        coords = f"{ref_name}:{ref_start}-{ref_end}"
        # Check if multimapped
        if aln.is_secondary:
            multimapped = 1
        else:
            multimapped = 0
    return [
        slice_name,
        parent_read,
        pe,
        slice_number,
        uid,
        mapped,
        multimapped,
        ref_name,
        ref_start,
        ref_end,
        coords,
    ]


@get_timing(task_name="processing BAM file")
def parse_bam(bam):
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
     bam: File name of bam file to process.

    Returns:
     pd.Dataframe: DataFrame with the columns listed above.

    """

    df_bam = pd.DataFrame(
        [
            parse_alignment(aln)
            for aln in pysam.AlignmentFile(bam, "rb").fetch(until_eof=True)
        ],
        columns=[
            "slice_name",
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

    df_bam["parent_id"] = hash_column(df_bam["parent_read"])
    df_bam["parent_id"] = df_bam["parent_id"].astype(str)
    df_bam["chrom"] = df_bam["chrom"].astype("category")
    df_bam["pe"] = df_bam["pe"].astype("category")
    df_bam["slice"] = df_bam["slice"].astype(int)
    df_bam["uid"] = df_bam["uid"].astype(int)
    df_bam["start"] = df_bam["start"].replace("", 0).astype(int)
    df_bam["end"] = df_bam["end"].replace("", 0).astype(int)
    df_bam.set_index(["slice_name", "chrom", "start"], inplace=True)
    return df_bam


class CCHDF5ReaderProcess(multiprocessing.Process):
    def __init__(
        self,
        path: os.PathLike,
        key: str,
        inq: multiprocessing.Queue,
        outq: multiprocessing.Queue,
    ):

        # Reading vars
        self.path = path
        self.key = key

        # Multiprocessing vars
        self.inq = inq
        self.outq = outq
        self.block_period = 0.01

        super(CCHDF5ReaderProcess, self).__init__()

    def run(self):

        with pd.HDFStore(self.path, "r") as store:
            while True:
                try:
                    viewpoint_to_find = self.inq.get(
                        block=True, timeout=self.block_period
                    )

                    if viewpoint_to_find is None:
                        break

                    df = store.select(
                        self.key,
                        where="viewpoint in viewpoint_to_find",
                        columns=["parent_id", "restriction_fragment"],
                    )

                    if not df.empty:
                        self.outq.put((viewpoint_to_find, df))

                except queue.Empty:
                    pass


class FragmentCountingProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        outq: multiprocessing.Queue,
    ):

        # Multiprocessing vars
        self.inq = inq
        self.outq = outq
        self.block_period = 0.01

        super(FragmentCountingProcess, self).__init__()

    def run(self):

        while True:
            try:
                vp, df = self.inq.get(block=True, timeout=self.block_period)

                if df is None:
                    break

                counts = (
                    df.groupby(["parent_id"])
                    .apply(get_fragment_combinations)
                    .reset_index(drop=True)
                    .explode()
                    .to_frame("combinations")
                    .assign(
                        bin1_id=lambda df: df["combinations"].map(lambda c: c[0]),
                        bin2_id=lambda df: df["combinations"].map(lambda c: c[1]),
                    )
                    .groupby(["bin1_id", "bin2_id"])
                    .size()
                    .reset_index()
                    .rename(columns={0: "count"})
                )

                self.outq.put((vp, counts))

            except queue.Empty:
                pass


class CCHDF5WriterProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        output_path: os.PathLike,
        output_key: str = "/",
        output_format: str = "cooler",
        single_file: bool = True,
        restriction_fragment_map: os.PathLike = None,
        viewpoint_path: os.PathLike = None,
    ):

        # Writing vars
        self.path = output_path
        self.key = output_key
        self.output_format = output_format
        self.single_file = single_file
        self.restriction_fragment_map = restriction_fragment_map
        self.viewpoint_path = viewpoint_path

        # Multiprocessing vars
        self.inq = inq
        self.block_period = 0.01

        super(CCHDF5WriterProcess, self).__init__()

    def run(self):

        while True:
            try:
                vp, df = self.inq.get(block=True, timeout=self.block_period)

                if df is None:
                    break

                if self.output_format == "tsv" and self.single_file:
                    df.to_csv(self.path, sep="\t", mode="a")

                elif self.output_format == "hdf5":
                    df.to_hdf(self.path, self.key, format="table", mode="a")

                elif self.output_format == "cooler":

                    if self.restriction_fragment_map and self.viewpoint_path:

                        from capcruncher.tools.storage import create_cooler_cc

                        create_cooler_cc(
                            output_prefix=self.path,
                            pixels=df,
                            bins=self.restriction_fragment_map,
                            viewpoint_name=vp,
                            viewpoint_path=self.viewpoint_path,
                        )
                    else:
                        raise ValueError(
                            "Restriction fragment map or path to viewpoints not supplied"
                        )

            except queue.Empty:
                pass
