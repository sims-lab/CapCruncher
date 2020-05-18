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
p.add_argument('-l', '--logfile', help='filename for logfile', default=sys.stdout)
p.add_argument('--buffer', help='Number of reads to process before writing output', default=10000, type=int)
args = p.parse_args()

def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return open(fn, 'w')
    else:
        return fn

def read_paired_fastq(fq1, fq2, read_counter, outq):
    '''Reads R1 and R2 fastq files and places paired reads into a queue'''
    counter = 0
    for r1, r2 in zip(fq1, fq2):
        counter += 1
        outq.put((r1, r2))
        
        if counter % 10000 == 0:
            print(f'Processed {counter} reads')
    
    outq.put((None, None))
    read_counter.value = counter

def remove_read_duplicates(inq, outq, reads_removed):
    '''Hashes the combined read1/read2 sequence and discards duplicates'''
    seen = set()
    removed_counter = 0
    r1, r2 = inq.get() # Gets reads from the input queue
    while r1:
        read_pair = hash(r1.sequence + r2.sequence)
        
        if read_pair not in seen:
            seen.add(read_pair)
            outq.put((r1,r2)) # Puts deduplicated reads in output queue
        else:
            removed_counter += 1
              
        r1, r2 = inq.get()
        
    outq.put((None, None))
    reads_removed.value = removed_counter


def write_to_fastq(inq):
    '''Writes all deduplicated read1/read2 to the appropriate file'''
    counter = 1
    r1_buffer, r2_buffer = [], []
    with xopen(filename=args.out1, mode='wb', compresslevel=args.compression_level, threads=4) as f1,\
         xopen(filename=args.out2, mode='wb', compresslevel=args.compression_level, threads=4) as f2:
        
        r1, r2 = inq.get()
        while r1:
            r1_buffer.append(str(r1))
            r2_buffer.append(str(r2))
            
            if counter % args.buffer == 0:
                f1.write('\n'.join(r1_buffer).encode())
                f2.write('\n'.join(r2_buffer).encode())
                r1_buffer, r2_buffer = [], []
            
            r1, r2 = inq.get()
            counter += 1
        
        # If the queue is empty write the rest
        f1.write('\n'.join(r1_buffer).encode())
        f2.write('\n'.join(r2_buffer).encode())
        
        
            
def main():
    
    # Set up multiprocessing variables 
    q1 = mp.SimpleQueue() # reads are placed into this queue for deduplication
    q2 = mp.SimpleQueue() # deduplicated reads are placed into the queue for writing
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
