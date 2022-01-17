import glob
import logging
import multiprocessing
import os
import pathlib
import queue
import random
import string
import traceback
from typing import Literal, Union

import numpy as np
import pandas as pd
import pysam
import tqdm
from capcruncher.tools.count import (
    get_fragment_combinations,
    preprocess_reporters_for_counting,
)
from capcruncher.utils import get_timing, hash_column
from pysam import FastxFile
from xopen import xopen
import xxhash
from collections import namedtuple


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
            for (read_counter, read) in enumerate(zip(*self._input_files_pysam)):

                # print(f"read_counter: {read_counter}, read: {read}, read_buffer: {self.read_buffer}")
                buffer.append(read)
                if read_counter % self.read_buffer == 0 and not read_counter == 0:
                    self.outq.put(buffer.copy())
                    buffer.clear()
                    logging.info(f"{read_counter} reads parsed (batch)")
                    rc = read_counter
                else:
                    rc = read_counter

            self.outq.put(buffer)  # Deal with remainder
            self.outq.put_nowait(None)  # Poison pill to terminate queue
            logging.info(f"{rc} reads parsed (final)")

        except Exception as e:
            logging.info(f"Reader failed with exception: {e}")
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
    ):

        super(FastqWriterProcess, self).__init__()

        self.inq = inq
        self.output = output
        self.compression_level = compression_level
        self.file_handles = self._get_filehandles()
        self.name = "FastqWriter"

    def _get_filehandles(self):
        if isinstance(self.output, str):
            return [
                xopen(
                    self.output, "w", compresslevel=self.compression_level, threads=0
                ),
            ]
        elif isinstance(self.output, (list, tuple, pd.Series)):
            return [
                xopen(fn, "w", compresslevel=self.compression_level, threads=0)
                for fn in self.output
            ]

    def _inputs_match_number_of_handles(self, reads):
        if len(reads[0]) == len(self.output):
            return True

    def run(self):

        while True:
            try:
                reads = self.inq.get(block=True, timeout=0.01)
                if reads:
                    is_string_input = True if isinstance(reads, str) else False

                    if is_string_input:
                        for fh in self.file_handles:
                            fh.write(reads)

                    else:
                        reads_str = [
                            "\n".join([str(r) for r in read_glob])
                            for read_glob in zip(*reads)
                        ]

                        for fh, read_set in zip(self.file_handles, reads_str):
                            fh.write((read_set + "\n"))

                else:

                    for fh in self.file_handles:
                        fh.close()

                    break

            except queue.Empty:
                continue


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


def parse_alignment(aln) -> CCAlignment:
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
        pe=pe,
        slice=int(slice_number),
        uid=int(uid),
        mapped=mapped,
        multimapped=multimapped,
        chrom=ref_name,
        start=int(ref_start),
        end=int(ref_end),
        coordinates=coords,
    )


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

    # Load reads into dataframe
    df_bam = pd.DataFrame(
        [
            parse_alignment(aln)
            for aln in pysam.AlignmentFile(bam, "rb").fetch(until_eof=True)
        ],
    )

    # Perform dtype conversions
    df_bam["chrom"] = df_bam["chrom"].astype("category")
    pe_category = pd.CategoricalDtype(["flashed", "pe"])
    df_bam["pe"] = df_bam["pe"].astype(
        pe_category
    )  # Only the one type present so need to include both

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

    def _select_by_viewpoint(self, store, viewpoint):

        viewpoint = viewpoint
        df = store.select(
            self.key,
            where="viewpoint in viewpoint",
            columns=[
                "parent_id",
                "restriction_fragment",
                "viewpoint",
                "capture",
                "exclusion",
            ],
        )
        return df

    def run(self):

        with pd.HDFStore(self.path, "r") as store:
            while True:
                try:
                    viewpoint_to_find = self.inq.get(
                        block=True, timeout=self.block_period
                    )

                    if viewpoint_to_find is None:
                        break

                    df = self._select_by_viewpoint(store, viewpoint_to_find)
                    if not df.empty:
                        self.outq.put((viewpoint_to_find, df))

                except queue.Empty:
                    pass


class CCParquetReaderProcess(multiprocessing.Process):
    def __init__(
        self,
        path: os.PathLike,
        inq: multiprocessing.Queue,
        outq: multiprocessing.Queue,
        selection_mode: Literal["single", "batch", "partition"],
    ):

        # Reading vars
        self.path = path
        self.partitions = glob.glob(f"{path}/*.parquet") or [
            self.path,
        ]
        self.selection_mode = selection_mode

        # Multiprocessing vars
        self.inq = inq
        self.outq = outq
        self.block_period = 0.01

        super(CCParquetReaderProcess, self).__init__()

    def _select_by_viewpoint(self, viewpoint):

        df = pd.read_parquet(
            self.path,
            columns=[
                "restriction_fragment",
                "viewpoint",
                "capture",
                "exclusion",
            ],
            filters=[[("viewpoint", "==", viewpoint)]],
            engine="pyarrow",
        )
        return df

    def _select_by_viewpoint_batch(self, viewpoints):

        df = pd.read_parquet(
            self.path,
            columns=[
                "parent_id",
                "restriction_fragment",
                "viewpoint",
                "capture",
                "exclusion",
            ],
            filters=[("viewpoint", "in", viewpoints)],
            engine="pyarrow",
        )
        return df

    def _select_by_viewpoint_batch_by_partition(self, viewpoints):

        for part in self.partitions:
            df = pd.read_parquet(
                part,
                columns=[
                    "restriction_fragment",
                    "viewpoint",
                    "capture",
                    "exclusion",
                ],
                filters=[("viewpoint", "in", viewpoints)],
                engine="pyarrow",
            )
            yield df

    def run(self):

        while True:
            try:
                viewpoints_to_find = self.inq.get(block=True, timeout=self.block_period)

                if viewpoints_to_find is None:
                    break

                elif (
                    self.selection_mode == "single"
                ):  # Slower as need to read all partitions each time
                    assert isinstance(viewpoints_to_find, str)
                    df = self._select_by_viewpoint(viewpoints_to_find)
                    if not df.empty:
                        self.outq.put((viewpoints_to_find, df))

                elif (
                    self.selection_mode == "batch"
                ):  # Faster as very low overhead when increasing number of viewpoints
                    df = self._select_by_viewpoint_batch(viewpoints_to_find)
                    for vp, df_vp in df.groupby("viewpoint"):
                        if not df_vp.empty:
                            self.outq.put((vp, df_vp))

                elif (
                    self.selection_mode == "partition"
                ):  # Low memory counting. Good for data rich viewpoints
                    for df in self._select_by_viewpoint_batch_by_partition(
                        viewpoints_to_find
                    ):
                        for vp, df_vp in df.groupby("viewpoint"):
                            if not df_vp.empty:
                                self.outq.put((vp, df_vp))

            except queue.Empty:
                pass


class FragmentCountingProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        outq: multiprocessing.Queue,
        **preprocessing_kwargs,
    ):

        # Multiprocessing vars
        self.inq = inq
        self.outq = outq
        self.block_period = 0.01
        self.preprocessing_kwargs = preprocessing_kwargs

        super(FragmentCountingProcess, self).__init__()

    def run(self):

        while True:
            try:
                vp, df = self.inq.get(block=True, timeout=self.block_period)

                if df is None:
                    break

                df = preprocess_reporters_for_counting(df, **self.preprocessing_kwargs)

                counts = (
                    df.groupby(["parent_id"])
                    .apply(get_fragment_combinations)
                    .reset_index(drop=True)
                    .explode()
                    .to_frame("combinations")
                    .dropna(axis=0)
                    .assign(
                        bin1_id=lambda df: df["combinations"].map(lambda c: c[0]),
                        bin2_id=lambda df: df["combinations"].map(lambda c: c[1]),
                    )
                    .groupby(["bin1_id", "bin2_id"])
                    .size()
                    .reset_index()
                    .rename(columns={0: "count"})
                    .sort_values(["bin1_id", "bin2_id"])
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
        n_viewpoints: int = 0,
    ):

        # Writing vars
        self.path = output_path
        self.key = output_key
        self.output_format = output_format
        self.single_file = single_file
        self.restriction_fragment_map = restriction_fragment_map
        self.viewpoint_path = viewpoint_path
        self.n_viewpoints = n_viewpoints

        if self.output_format == "cooler":
            self.bins = pd.read_csv(
                restriction_fragment_map,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "name"],
            )

        # Multiprocessing vars
        self.inq = inq
        self.block_period = 0.01

        super(CCHDF5WriterProcess, self).__init__()

    def run(self):

        with tqdm.tqdm(total=self.n_viewpoints) as pbar:

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
                                bins=self.bins,
                                viewpoint_name=vp,
                                viewpoint_path=self.viewpoint_path,
                                ordered=True,
                                dupcheck=False,
                                triucheck=False,
                            )
                        else:
                            raise ValueError(
                                "Restriction fragment map or path to viewpoints not supplied"
                            )

                    pbar.update()

                except queue.Empty:
                    pass


class CCCountsWriterProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        output_format: str = "cooler",
        restriction_fragment_map: os.PathLike = None,
        viewpoint_path: os.PathLike = None,
        tmpdir: os.PathLike = ".",
    ):

        # Writing vars
        self.output_format = output_format
        self.restriction_fragment_map = restriction_fragment_map
        self.viewpoint_path = viewpoint_path
        self.tmpdir = tmpdir

        if self.output_format == "cooler":
            self.bins = pd.read_csv(
                restriction_fragment_map,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "name"],
            )

        # Multiprocessing vars
        self.inq = inq
        self.block_period = 0.01

        super(CCCountsWriterProcess, self).__init__()

    def run(self):

        while True:
            try:
                vp, df = self.inq.get(block=True, timeout=self.block_period)

                if df is None:
                    break
                else:
                    path = os.path.join(
                        self.tmpdir,
                        "".join(
                            random.choices(string.ascii_uppercase + string.digits, k=6)
                        ),
                    )

                    if self.output_format == "tsv" and self.single_file:
                        df.to_csv(path, sep="\t", mode="a")

                    elif self.output_format == "hdf5":
                        df.to_hdf(path, vp, format="table", mode="a")

                    elif self.output_format == "cooler":

                        if self.restriction_fragment_map and self.viewpoint_path:

                            from capcruncher.tools.storage import create_cooler_cc

                            create_cooler_cc(
                                output_prefix=path,
                                pixels=df,
                                bins=self.bins,
                                viewpoint_name=vp,
                                viewpoint_path=self.viewpoint_path,
                                ordered=True,
                                dupcheck=False,
                                triucheck=False,
                            )
                        else:
                            raise ValueError(
                                "Restriction fragment map or path to viewpoints not supplied"
                            )

            except queue.Empty:
                pass
