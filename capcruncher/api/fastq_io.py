from collections.abc import Callable, Sequence
import multiprocessing
import os
from pathlib import Path
from typing import Any, cast

from loguru import logger

from pysam import FastxFile
from xopen import xopen

type FastqFormatFunction = Callable[[object], object]


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
        input_files: Path | str | Sequence[Path | str],
        outq: multiprocessing.Queue,
        read_buffer: int = 100000,
    ) -> None:
        # Input variables
        self.input_files = self._normalise_input_files(input_files)
        self._multifile = len(self.input_files) > 1

        # Multiprocessing variables
        self.outq = outq

        # Reader variables
        self.read_buffer = read_buffer

        super().__init__()

    def _normalise_input_files(
        self, input_files: Path | str | Sequence[Path | str]
    ) -> list[str]:
        if isinstance(input_files, str | os.PathLike):
            return [os.fspath(input_files)]
        return [os.fspath(input_file) for input_file in input_files]

    def run(self) -> None:
        """Performs reading and chunking of fastq file(s)."""

        if self._multifile:
            input_files_pysam = [FastxFile(f) for f in self.input_files]
        else:
            input_files_pysam = [
                FastxFile(self.input_files[0]),
            ]

        try:
            buffer = []
            rc = 0
            for read_counter, read in enumerate(zip(*input_files_pysam)):
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
            self.outq.put("END")  # Poison pill to terminate queue
            logger.info(f"{rc} reads parsed (final)")

        except Exception as e:
            logger.info(f"Reader failed with exception: {e}")
            raise

        finally:
            for fh in input_files_pysam:
                fh.close()


class FastqReadFormatterProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        outq: multiprocessing.Queue,
        formatting: Sequence[FastqFormatFunction] | None = None,
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

        super().__init__()

    def _format_as_str(self, reads: Sequence[Sequence[object]]) -> list[str]:
        # [(r1, r2), (r1, r2)] -> [r1 combined string, r2 combined string]
        return ["\n".join([str(rn) for rn in r]) for r in zip(*reads)]

    def run(self) -> None:
        try:
            reads = self.inq.get()

            while not reads == "END":
                for formatting_to_apply in self.formatting:
                    reads = formatting_to_apply(
                        cast(Sequence[Sequence[object]], reads)
                    )

                self.outq.put(reads)
                reads = self.inq.get()

            self.outq.put("END")

        except Exception:
            logger.exception("Formatter worker failed")
            self.outq.put("END")


class FastqWriterSplitterProcess(multiprocessing.Process):
    def __init__(
        self,
        inq: multiprocessing.Queue,
        output_prefix: Path | str,
        paired_output: bool = False,
        gzip: bool = False,
        compression_level: int = 3,
        compression_threads: int = 8,
        n_subprocesses: int = 1,
        n_workers_terminated: int = 0,
        n_files_written: int = 0,
    ) -> None:
        self.inq = inq
        self.output_prefix = os.fspath(output_prefix)
        self.paired_output = paired_output

        self.gzip = gzip
        self.compression_level = compression_level
        self.compression_threads = compression_threads

        self.n_subprocesses = n_subprocesses
        self.n_workers_terminated = n_workers_terminated
        self.n_files_written = n_files_written

        super().__init__()

    def _get_file_handles(self) -> list[Any]:
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

    def run(self) -> None:
        try:
            reads = self.inq.get()
            is_string_input = True if isinstance(reads[0], str) else False

            while self.n_workers_terminated < self.n_subprocesses:
                if reads == "END":
                    self.n_workers_terminated += 1
                    continue

                elif is_string_input:
                    for fh, read in zip(self._get_file_handles(), reads):
                        fh.write(read + "\n")
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
            logger.exception("Writer worker failed")
