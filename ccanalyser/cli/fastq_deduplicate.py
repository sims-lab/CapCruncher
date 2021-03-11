#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""
import os
import re
from collections import Counter
from multiprocessing import SimpleQueue
from typing import List, Tuple, Union

import click
import numpy as np
import pandas as pd
import ujson
from ccanalyser.cli import cli
from ccanalyser.tools.deduplicate import (ReadDeduplicationParserProcess,
                                          ReadDuplicateRemovalProcess)
from ccanalyser.tools.io import FastqReaderProcess, FastqWriterProcess
from ccanalyser.tools.statistics import DeduplicationStatistics
from ccanalyser.utils import NaturalOrderGroup, invert_dict, load_json
from xopen import xopen


@cli.group(cls=NaturalOrderGroup)
@click.pass_context
def fastq_deduplicate(ctx):
    """Identifies PCR duplicate fragments from Fastq files"""


@fastq_deduplicate.command()
@click.argument("input_files", nargs=-1)
@click.option(
    "-o",
    "--output",
    help="File to store hashed sequence identifiers",
    default="out.json",
)
@click.option(
    "--read_buffer",
    help="Number of reads to process before writing to file",
    default=1e5,
    type=click.INT,
)
def parse(input_files: Tuple, output: os.PathLike = "out.json", read_buffer: int = 1e5):
    """
    Parses fastq file(s) into easy to deduplicate format.

    Generates a hashed (using xxhash) dictionary in json format for deduplication using identify. 

    Args:
     input_files (Tuple): One or more fastq files to process
     output (os.PathLike, optional): Output for parsed read identifiers and sequences. Defaults to "out.json".
     read_buffer (int, optional): Number of reads to process before outputting to file. Defaults to 1e5.
    """    


    # Set up multiprocessing variables
    inputq = SimpleQueue()  # Reads are placed into this queue for deduplication
    writeq = SimpleQueue()  # Deduplicated reads are placed into the queue for writing

    reader = FastqReaderProcess(
        input_files=input_files,
        outq=inputq,
        n_subprocesses=1,
        read_buffer=read_buffer,
    )

    parser = ReadDeduplicationParserProcess(
        inq=inputq, outq=writeq, save_hashed_dict_path=output
    )

    processes = [reader, parser]

    for proc in processes:
        proc.start()

    for proc in processes:
        proc.join()
        proc.terminate()


@fastq_deduplicate.command()
@click.argument(
    "input_files",
    nargs=-1,
)
@click.option(
    "-o", "--output", help="Output file", default="duplicates.json", required=True
)
def identify(input_files: Tuple, output: os.PathLike = "duplicates.json"):
    """Identifies fragments with duplicated sequences.

       Takes the input of the parse subcommand and deduplicates the hashed sequences 
       using a python dictionary. Duplicate read ids (hashed) are output  as a .json file.

    Args:
     input_files (Tuple): Paths to json files containing dictionaries with hashed read ids as the keys
                          and hashed sequences as the values.
     output (os.PathLike, optional): Duplicate read ids identified. Defaults to "duplicates.json".
    """    


    dedup_sequences = dict()
    read_ids = set()

    np.random.shuffle(np.array(input_files))
    for fn in input_files:
        d = load_json(fn)  # {READ_NAME_HASH: SEQUENCE_HASH}
        read_ids.update(d)
        dedup_sequences.update(invert_dict(d))  # {SEQUENCE_HASH: READ_NAME_HASH}

    duplicated_ids = read_ids - set(dedup_sequences.values())
    del read_ids
    del dedup_sequences

    with xopen(output, "w") as w:
        duplicated_ids_dict = dict.fromkeys(duplicated_ids)
        ujson.dump(duplicated_ids_dict, w)


@fastq_deduplicate.command()
@click.argument("input_files", nargs=-1)
@click.option(
    "-o",
    "--output_prefix",
    help="Output prefix for deduplicated fastq file(s)",
    default="",
)
@click.option(
    "-d",
    "--duplicated_ids",
    help="Path to duplicate ids, identified by the identify subcommand",
)
@click.option(
    "--read_buffer",
    help="Number of reads to process before writing to file",
    default=1e5,
    type=click.INT,
)
@click.option(
    "--gzip/--no-gzip",
    help="Determines if files are gziped or not",
    default=False
)
@click.option(
    "--compression_level",
    help="Level of compression for output files",
    default=5,
    type=click.INT,
)
@click.option("--sample_name", help="Name of sample e.g. DOX_treated_1", default='sampleX')
@click.option("--stats_prefix", help="Output prefix for stats file", default='stats')
def remove(
    input_files: Tuple,
    duplicated_ids: os.PathLike,
    read_buffer: int = 1e5,
    output_prefix: os.PathLike = "",
    gzip: bool = False, 
    compression_level: int = 5,
    sample_name: str = "",
    stats_prefix: os.PathLike="",
):
    """
    Removes fragments with duplicated sequences from fastq files.
    
    Reads fastq files and removes all duplicate entries provided by the duplicated ids file.
    Uses the output of the identify subcommand to identify and remove all sequence duplicates.

    Args:
     input_files (Tuple): Input fastq files (in the same order as used for the parse command).
     duplicated_ids (os.PathLike): Duplicated read ids from identify command (hashed and in json format).
     read_buffer (int, optional): Number of reads to process before writing to file. Defaults to 1e5.
     output_prefix (os.PathLike, optional): Deduplicated fastq output prefix. Defaults to "".
     gzip (bool, optional): Determines if output is gzip compressed using pigz. Defaults to False.
     compression_level (int, optional): Level of compression if required (1-9). Defaults to 5.
     sample_name (str, optional): Name of sample processed e.g. DOX-treated_1. Defaults to "".
     stats_prefix (os.PathLike, optional): Output prefix for statistics. Defaults to "".
    
    """

    duplicated_ids = set(load_json(duplicated_ids))
    inputq = SimpleQueue()  # Reads are placed into this queue for deduplication
    writeq = SimpleQueue()  # Deduplicated reads are placed into the queue for writing
    statq = SimpleQueue()  # Statistics are sent on this queue for processing

    output_files = [
        f"{output_prefix}_{ii+1}.fastq{'.gz' if gzip else ''}" for ii in range(len(input_files))
    ]

    deduplicator = [
        ReadDuplicateRemovalProcess(
            inq=inputq, outq=writeq, duplicated_ids=duplicated_ids, statq=statq
        )
        for _ in range(1) # Placeholder to enable multiple digestion processes at a later date
    ]

    del duplicated_ids # Reduces memory usage before starting (likely by forking) a new process

    reader = FastqReaderProcess(
        input_files=input_files,
        outq=inputq,
        read_buffer=read_buffer,
        n_subprocesses=1,
    )

    writer = FastqWriterProcess(
        inq=writeq,
        output=output_files,
        compression_level=compression_level,
    )

    reader.start()
    writer.start()
    for dedup in deduplicator:
        dedup.start()

    processes = [writer, reader, *deduplicator]

    for proc in processes:
        proc.join()
        proc.terminate()

    # Handle statistics
    stats_aggregator = Counter()
    stats = statq.get()

    while not stats == "END":
        stats_aggregator.update(stats)
        stats = statq.get()

    deduplication_stats = DeduplicationStatistics(
        sample=sample_name, **stats_aggregator
    )

    deduplication_stats.df.to_csv(f"{stats_prefix}.deduplication.csv", index=False)
