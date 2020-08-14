import multiprocessing as mp
import os
import re
import sys
from collections import Counter
import time

import pandas as pd
import numpy as np
import pysam
from pysam import FastxFile
from xopen import xopen
from .digest_genome import get_re_site


class DigestedRead:
    def __init__(
        self,
        read,
        cutsite,
        min_slice_length=18,
        slice_number_offset=0,
        allow_undigested=False,
        read_type="flashed",
    ):

        self.read = read
        self.min_slice_length = min_slice_length
        self.slice_number_offset = slice_number_offset
        self.allow_undigested = allow_undigested
        self.read_type = read_type

        self.recognition_seq = cutsite.upper()
        self.recognition_len = len(cutsite)
        self.recognition_re = re.compile(self.recognition_seq)

        self.slice_indexes = self.get_recognition_site_indexes()
        self.slices_total_counter = len(self.slice_indexes) - 1
        self.slices_valid_counter = 0
        self.has_slices = self.slices_total_counter > 1

    def get_recognition_site_indexes(self):
        indexes = [
            re_site.start()
            for re_site in self.recognition_re.finditer(self.read.sequence.upper())
        ]

        indexes.insert(0, 0)
        indexes.append(len(self.read.sequence))

        return indexes

    @property
    def slices(self):

        indexes = self.slice_indexes
        slice_no = 0 + self.slice_number_offset
        slices_list = []

        if self.has_slices or self.allow_undigested:

            for ii, (slice_start, slice_end) in enumerate(zip(indexes, indexes[1:])):

                if ii > 0:
                    slice_start += self.recognition_len

                if self.is_valid_slice(slice_start, slice_end):
                    slices_list.append(
                        self.prepare_slice(slice_start, slice_end, slice_no)
                    )

                    self.slices_valid_counter += 1
                    slice_no += 1

        return slices_list

    def prepare_slice(self, start, end, slice_no):
        return "\n".join(
            [
                f"@{self.read.name}{self.read.comment}|{self.read_type}|{slice_no}",
                self.read.sequence[start:end],
                "+",
                self.read.quality[start:end],
            ]
        )

    def is_valid_slice(self, start, end):
        if end - start >= self.min_slice_length:
            return True

    def __str__(self):
        slices = self.slices
        return ("\n".join(slices) + "\n") if slices else ""

def read_fastq(fq, outq, n_workers=1, buffer=10000):

    """Reads fastq file and places reads into a queue"""
    r_buffer = []
    for rc, read in enumerate(fq):
        r_buffer.append(read)

        if rc % buffer == 0 and not rc == 0:
            outq.put(r_buffer)
            r_buffer = []
            print(f"Processed {rc} reads")

    outq.put(r_buffer)  # Add the reads that don't fit the buffer size

    for _ in range(n_workers):
        outq.put("TER")


def read_paired_fastqs(fq1, fq2, outq, n_workers=1, buffer=10000):
    """Reads R1 and R2 fastq files and places paired reads into a queue"""
    r_buffer = []
    for rc, (r1, r2) in enumerate(zip(fq1, fq2)):
        r_buffer.append((r1, r2))

        if rc % buffer == 0 and not rc == 0:
            outq.put(r_buffer)
            r_buffer = []
            print(f"Processed {rc} reads")

    outq.put(r_buffer)
    for _ in range(n_workers):
        outq.put("TER")  # Places a terminator in the queue for every spawned worker


def digest_read_flashed(inq, outq, statq, **kwargs):

    read_buffer = []
    stat_buffer = []
    reads = inq.get()

    while not reads == "TER":
        for read in reads:
            sliced_read = DigestedRead(read, **kwargs)
            sliced_read_str = str(sliced_read)
            if sliced_read_str:  # Only append if valid slices present
                read_buffer.append(sliced_read_str)

            stat_buffer.append(
                (sliced_read.slices_total_counter, sliced_read.slices_valid_counter,)
            )

        outq.put(read_buffer)  # Add list of digested reads to the queue
        statq.put(stat_buffer)
        read_buffer = []  # Clear buffer
        stat_buffer = []  # Clear buffer
        reads = inq.get()  # Get new list of reads

    outq.put("TER")
    statq.put("TER")


def digest_read_unflashed(inq, outq, statq, **kwargs):

    read_buffer = []
    stat_buffer = []
    reads = inq.get()

    while not reads == "TER":  # Checks to see if the queue has been terminated
        for r1, r2 in reads:
            sliced_read_1 = DigestedRead(r1, **kwargs)  # Digest read 1
            kwargs["slice_number_offset"] = sliced_read_1.slices_valid_counter  # Update slice offset

            sliced_read_2 = DigestedRead(r2, **kwargs)  # Digest read 2

            kwargs["slice_number_offset"] = 0  # Reset slice offset

            s1, s2 = str(sliced_read_1), str(sliced_read_2)

            if s1 and s2:  # Only store of both reads have valid slices
                read_buffer.append(f"{s1}{s2}")

            stat_buffer.append(
                (
                    sliced_read_1.slices_total_counter,
                    sliced_read_1.slices_valid_counter,
                    sliced_read_2.slices_total_counter,
                    sliced_read_2.slices_valid_counter,
                )
            )

        outq.put(read_buffer)
        statq.put(stat_buffer)
        read_buffer = []
        stat_buffer = []
        reads = inq.get()

    outq.put("TER")
    statq.put("TER")


def write_to_fastq(inq, output_file="out.fq.gz", compression_level=5, n_workers=1):

    """Writes all digested reads to the appropriate file"""
    with xopen(filename=output_file, mode="wb", compresslevel=compression_level) as f:

        ter_counter = 0
        while ter_counter < n_workers:

            reads = inq.get()
            if not reads == "TER":
                f.write("".join(reads).encode())
            else:
                ter_counter += 1


def collate_stats_flashed(inq, stats_file="out.log", n_workers=1):

    total_counter = Counter()
    valid_counter = Counter()
    ter_counter = 0

    while ter_counter < n_workers:

        counts = inq.get()
        if not counts == "TER":
            total, valid = zip(*counts)
            total_counter = total_counter + Counter(total)
            valid_counter = valid_counter + Counter(valid)
        else:
            ter_counter += 1

    dframes = []
    for key, values in zip(["total", "valid"], [total_counter, valid_counter]):

        dframes.append(
            pd.DataFrame.from_dict(values, orient="index")
            .reset_index()
            .rename(columns={"index": "bin", 0: "frequency"})
            .assign(stat=key)
        )

    df_stats = pd.concat(dframes)
    df_stats["read_type"] = "flashed"
    df_stats["bin"] = df_stats["bin"].astype(int)

    (df_stats.sort_values("bin").to_csv(stats_file, index=False))


def collate_stats_unflashed(inq, stats_file="out.log", n_workers=1):
    counters = {
        "r1": {"total": Counter(), "valid": Counter()},
        "r2": {"total": Counter(), "valid": Counter()},
    }
    ter_counter = 0

    while ter_counter < n_workers:
        counts = inq.get()
        if not counts == "TER":
            r1_total, r1_valid, r2_total, r2_valid = zip(*counts)
            counters["r1"]["total"] = counters["r1"]["total"] + Counter(r1_total)
            counters["r1"]["valid"] = counters["r1"]["valid"] + Counter(r1_valid)
            counters["r2"]["total"] = counters["r2"]["total"] + Counter(r2_total)
            counters["r2"]["valid"] = counters["r2"]["valid"] + Counter(r2_valid)
        else:
            ter_counter += 1

    dframes = []
    for read in counters:
        for key, values in counters[read].items():

            dframes.append(
                pd.DataFrame.from_dict(values, orient="index")
                .reset_index()
                .rename(columns={"index": "bin", 0: "frequency"})
                .assign(stat=key, read_type=read)
            )

    df_stats = pd.concat(dframes)
    df_stats["bin"] = df_stats["bin"].astype(int)

    (df_stats.sort_values("bin").to_csv(stats_file, index=False))


def main(
    subcommand,
    input_fastq=None,
    fq1=None,
    fq2=None,
    cut_sequence=None,
    restriction_enzyme=None,
    keep_cutsite=False,
    output_file="out.fastq.gz",
    minimum_slice_length=18,
    stats_file=None,
    compression_level=5,
    n_digestion_processes=1,
    buffer=10000,
):

    # Set up multiprocessing variables
    inputq = mp.SimpleQueue()  # reads are placed into this queue for processing
    writeq = mp.SimpleQueue()  # digested reads are placed into the queue for writing
    statq = mp.SimpleQueue()  # stats queue

    # Variables
    cut_site = get_re_site(cut_sequence, restriction_enzyme)

    if subcommand == "flashed":  # Checks the subsubcommand to see in which mode to run

        fq = FastxFile(input_fastq)

        # Writer process
        writer = mp.Process(
            target=write_to_fastq,
            args=(writeq,),
            kwargs={
                "output_file": output_file,
                "n_workers": n_digestion_processes,
                "compression_level": compression_level,
            },
        )
        writer.start()

        # Define processes
        processes = [
            mp.Process(
                target=read_fastq,
                args=(fq, inputq),
                kwargs={"n_workers": n_digestion_processes, "buffer": buffer},
            ),
            mp.Process(
                target=collate_stats_flashed,
                args=(statq,),
                kwargs={"stats_file": stats_file, "n_workers": n_digestion_processes},
            ),
        ]

        processes_repeated = [
            mp.Process(
                target=digest_read_flashed,
                args=(inputq, writeq, statq),
                kwargs={
                    "cutsite": cut_site,
                    "min_slice_length": minimum_slice_length,
                    'read_type': subcommand
                },
            )
            for i in range(n_digestion_processes)
        ]

        processes = processes + processes_repeated

    elif subcommand == "unflashed":

        fq1, fq2 = FastxFile(fq1), FastxFile(fq2)

        # Writer process
        writer = mp.Process(
            target=write_to_fastq,
            args=(writeq,),
            kwargs={"output_file": output_file, "n_workers": n_digestion_processes},
        )
        writer.start()

        # Define processes
        processes = [
            mp.Process(
                target=read_paired_fastqs,
                args=(fq1, fq2, inputq),
                kwargs={"n_workers": n_digestion_processes, "buffer": buffer},
            ),
            mp.Process(
                target=collate_stats_unflashed,
                args=(statq,),
                kwargs={"stats_file": stats_file, "n_workers": n_digestion_processes},
            ),
        ]

        processes_repeated = [
            mp.Process(
                target=digest_read_unflashed,
                args=(inputq, writeq, statq),
                kwargs={
                    "cutsite": cut_site,
                    "min_slice_length": minimum_slice_length,
                    "slice_number_offset": 0,
                    'allow_undigested': True,
                    'read_type': subcommand
                },
            )
            for i in range(n_digestion_processes)
        ]

        processes = processes + processes_repeated

    # Start all processes
    for proc in processes:
        proc.start()

    # Join processes (wait for processes to finish before the main process)
    writer.join()

    for proc in processes:
        proc.join(10)
        proc.terminate()
