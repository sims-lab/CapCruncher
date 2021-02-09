#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:45:09 2020

@author: asmith

Script splits a fastq into specified chunks
"""

from multiprocessing import Manager, SimpleQueue, Pipe
import click
from ccanalyser.cli import cli
from ccanalyser.tools.io import FastqReaderProcess, FastqWriterSplitterProcess, FastqReadFormatterProcess

@cli.command()
@click.argument("input_files", nargs=-1)
@click.option(
    "-o",
    "--output_prefix",
    help="Output prefix for deduplicated fastq file(s)",
    default="",
)
@click.option(
    "--compression_level",
    help="Level of compression for output files",
    default=5,
    type=click.INT,
)
@click.option(
    "-n",
    "--n_reads",
    help="Number of reads per fastq file",
    default=1e6,
    type=click.INT,
)
def fastq_split(input_files, output_prefix, compression_level=5, n_reads=1000000):

    '''Splits fastq file(s) into equal chunks'''

    readq = SimpleQueue()
    writeq = SimpleQueue()

    paired = True if len(input_files) > 1 else False

    reader = FastqReaderProcess(
        input_files=input_files, outq=readq, read_buffer=n_reads, n_subprocesses=1,
    )

    formatter = [FastqReadFormatterProcess(inq=readq, outq=writeq) 
                 for _ in range(1)]
    
    
    writer = FastqWriterSplitterProcess(
        inq=writeq, output_prefix=output_prefix, paired_output=paired, n_subprocesses=1)

    
    processes = [writer, reader, *formatter ]

    for proc in processes:
        proc.start()
    
    for proc in processes:
        proc.join()
        proc.terminate()
