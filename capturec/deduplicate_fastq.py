#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""
import os
import argparse
import sys
import multiprocessing as mp
from pysam import FastxFile
from xopen import xopen

# Make sure the script can find capturec scripts
SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(PACKAGE_DIR)

from bin.capturec_commandline import get_args_deduplicate_fastq as get_args


def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return open(fn, 'w')
    else:
        return fn

def read_paired_fastq(fq1,
                      fq2,
                      read_counter,
                      outq,
                      n_reads_buffer=10000):
    '''Reads R1 and R2 fastq files and places paired reads into a queue'''

    counter = 0
    buffer = []
    for r1, r2 in zip(fq1, fq2):
        counter += 1
        buffer.append((r1, r2))

        if counter % n_reads_buffer == 0:
            print(f'Processed {counter} reads')
            outq.put(buffer)
            buffer = []

    outq.put(buffer) # Add any reads that do not fit the batch size
    outq.put(None) # Used to terminate the queue
    read_counter.value = counter # Stores the number of processed reads.

def remove_read_duplicates(inq,
                           outq,
                           reads_removed,
                           n_reads_buffer=10000):
    '''Hashes the combined read1/read2 sequence and discards duplicates'''
    seen = set()
    removed_counter = 0
    buffer = []

    reads = inq.get() # Get list of reads from the input queue
    while reads:
        for read_counter, (r1, r2) in enumerate(reads):
            read_pair = hash(r1.sequence + r2.sequence)

            if read_pair not in seen:
                seen.add(read_pair)
                buffer.append( [str(r1), str(r2)] )
            else:
                removed_counter += 1

            if read_counter % n_reads_buffer == 0:
                outq.put(buffer)
                buffer = []


        reads = inq.get()

    outq.put(buffer) # Add any reads that do not fit the batch size
    outq.put(None) # Terminate the queue
    reads_removed.value = removed_counter # Store the number of removed reads


def write_to_fastq(inq,
                   out1,
                   out2,
                   compression_level=5):
    '''Writes all deduplicated read1/read2 to the appropriate file'''
    with xopen(filename=out1, mode='wb', compresslevel=compression_level) as f1,\
         xopen(filename=out2, mode='wb', compresslevel=compression_level) as f2:

        reads_paired = inq.get()
        while reads_paired:
            r1, r2 = zip(*reads_paired)
            f1.write(('\n'.join(r1) + '\n').encode())
            f2.write(('\n'.join(r2) + '\n').encode())
            reads_paired = inq.get()


def main(fq1,
         fq2,
         out1='out_1.fastq.gz',
         out2='out_2.fastq.gz',
         stats_file='stats_out.log',
         read_buffer=10000,
         compression_level=5):

    # Set up multiprocessing variables
    inputq = mp.SimpleQueue() # reads are placed into this queue for deduplication
    writeq = mp.SimpleQueue() # deduplicated reads are placed into the queue for writing
    manager = mp.Manager()
    r_counter = manager.Value('i', 0) # counts the number of readpairs
    d_counter = manager.Value('i', 0) # counts the number of duplicated readpairs

    # Fastq files
    fq1 = FastxFile(fq1)
    fq2 = FastxFile(fq2)

    # Writer process
    writer = mp.Process(target=write_to_fastq,
                        args=(writeq, out1, out2),
                        kwargs={'compression_level': compression_level})
    writer.start()

    # Define processes
    processes = [mp.Process(target=read_paired_fastq,
                            args=(fq1, fq2, r_counter, inputq),
                            kwargs={'n_reads_buffer': read_buffer}),
                 mp.Process(target=remove_read_duplicates,
                            args=(inputq, writeq, d_counter),
                            kwargs={'n_reads_buffer': read_buffer}),
                ]

    # Start processes
    for proc in processes:
        proc.start()

    # Join processes (wait for processes to finish before the main process)
    writer.join()

    # Terminate all processes
    for proc in processes:
        proc.terminate()

    # Write to log file
    with open_logfile(stats_file) as logfile:
            logfile.write(f'Read_pairs_processed\t{r_counter.value}\n')
            logfile.write(f'Read_pairs_unique\t{r_counter.value - d_counter.value}\n')
            logfile.write(f'Read_pairs_removed\t{d_counter.value}\n')


if __name__ == '__main__':
    args = get_args().parse_args()
    main(**vars(args))
