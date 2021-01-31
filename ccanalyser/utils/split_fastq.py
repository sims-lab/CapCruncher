#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:45:09 2020

@author: asmith

Script splits a fastq into specified chunks
"""

from multiprocessing import Manager, SimpleQueue, Pipe
from ccanalyser.utils.io import FastqReaderProcess, FastqWriterSplitterProcess, FastqReadFormatterProcess

def main(input_files, output_prefix, compression_level=5, n_reads=1000000, n_subprocesses=1):

    readq = SimpleQueue()
    writeq = SimpleQueue()

    paired = True if len(input_files) > 1 else False

    reader = FastqReaderProcess(
        input_files=input_files, outq=readq, read_buffer=n_reads, n_subprocesses=n_subprocesses,
    )

    formatter = [FastqReadFormatterProcess(inq=readq, outq=writeq) 
                 for _ in range(n_subprocesses)]
    
    
    writer = FastqWriterSplitterProcess(
        inq=writeq, output_prefix=output_prefix, paired_output=paired, n_subprocesses=n_subprocesses)

    
    processes = [writer, reader, *formatter ]

    for proc in processes:
        proc.start()
    
    for proc in processes:
        proc.join()
        proc.terminate()
