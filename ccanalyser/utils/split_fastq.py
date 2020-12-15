#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:45:09 2020

@author: asmith

Script splits a fastq into specified chunks
"""

from multiprocessing import SimpleQueue
from ccanalyser.utils.io import FastqReaderProcess, FastqWriterSplitterProcess

def main(input_files, output_prefix, compression_level=5, n_reads=1000000):

    queue = SimpleQueue()
    paired = True if len(input_files) > 1 else False

    reader = FastqReaderProcess(
        input_files=input_files, outq=queue, read_buffer=n_reads
    )

    writer = FastqWriterSplitterProcess(
        inq=queue, output_prefix=output_prefix, paired_output=paired,
    )

    
    processes = [reader, writer]

    for proc in processes:
        proc.start()
    
    for proc in processes:
        proc.join()
        proc.terminate()
