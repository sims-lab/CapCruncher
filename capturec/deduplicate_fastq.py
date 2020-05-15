#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""

import argparse
import sys
import multiprocessing as mp
from pysam import FastxFile
from xopen import xopen

p = argparse.ArgumentParser()
p.add_argument('-1', '--fq1', help='fastq file to parse containing read 1')
p.add_argument('-2', '--fq2', help='fastq file to parse containing read 2')
p.add_argument('--out1', help='fastq file to parse containing read 2', default='out1.fq.gz')
p.add_argument('--out2', help='fastq file to parse containing read 2', default='out2.fq.gz')
p.add_argument('-c', '--compression_level', help='Level of compression (1-9 with 9 being the highest)', type=int, default=5)
p.add_argument('-l', '--logfile', help='filename for logfile',
               default=sys.stdout)
args = p.parse_args()

def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return open(fn, 'w')
    else:
        return fn

def read_paired_fastq(fq1, fq2, read_counter, outq):
    '''Reads R1 and R2 fastq files and places paired reads into a queue'''
    for r1, r2 in zip(fq1, fq2):
        read_counter.value = read_counter.value + 1
        outq.put((r1, r2))
        
        if read_counter.value % 10000 == 0:
            print(f'Processed {read_counter.value} reads')
    
    outq.put((None, None))

def remove_read_duplicates(inq, outq, reads_removed):
    '''Hashes the combined read1/read2 sequence and discards duplicates'''
    seen = set()
    counter = 0
    r1, r2 = inq.get()
    while r1:
        read_pair = hash(r1.sequence + r2.sequence)
        
        if read_pair not in seen:
            seen.add(read_pair)
            outq.put((r1,r2))
        else:
            reads_removed.value = (reads_removed.value + 1)
              
        r1, r2 = inq.get()
        
    outq.put((None, None))


def write_to_fastq(inq):
    '''Writes all deduplicated read1/read2 to the appropriate file'''
    counter = 0
    with xopen(filename=args.out1, mode='wb', compresslevel=args.compression_level, threads=4) as f1,\
         xopen(filename=args.out2, mode='wb', compresslevel=args.compression_level, threads=4) as f2:
        
        r1, r2 = inq.get()
        while r1:
            f1.write((str(r1) + '\n' ).encode())
            f2.write((str(r2) + '\n').encode())
            
            r1, r2 = inq.get()
            
def main():
    
    # Set up multiprocessing variables 
    q1 = mp.Queue() # reads are placed into this queue for deduplication
    q2 = mp.Queue() # deduplicated reads are placed into the queue for writing
    manager = mp.Manager()
    r_counter = manager.Value('i', 0) # counts the number of readpairs
    d_counter = manager.Value('i', 0) # counts the number of duplicated readpairs

    # Fastq files
    fq1 = FastxFile(args.fq1)
    fq2 = FastxFile(args.fq2)

    # Writer process
    writer = mp.Process(target=write_to_fastq, args=(q2,))
    writer.start()

    # Define processes
    processes = [mp.Process(target=read_paired_fastq, args=(fq1, fq2, r_counter, q1)),
                 mp.Process(target=remove_read_duplicates, args=(q1, q2, d_counter)),
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
    with open_logfile(args.logfile) as logfile:
            logfile.write(f'Read_pairs_processed\t{r_counter.value}\n')
            logfile.write(f'Read_pairs_unique\t{r_counter.value - d_counter.value}\n')
            logfile.write(f'Read_pairs_removed\t{d_counter.value}\n')
         
    
if __name__ == '__main__':
    main()
