#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""

import argparse
import sys
from multiprocessing import Manager
from pysam import FastxFile
import gzip
from functools import partial
import itertools

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

def remove_duplicates(read_counter, read_encountered, reads_removed, reads):

    read_1, read_2 = reads
    read_pair = hash(read_1.sequence +  read_2.sequence)
    if not read_pair in read_encountered:
        read_encountered[read_pair] = 1
        return read_1, read_2
    else:
        reads_removed += 1

        
def main():

    fq1, fq2 = FastxFile(args.fq1), FastxFile(args.fq2)
    out_1 = gzip.open(args.out1, 'wb', compresslevel=args.compression_level)
    out_2 = gzip.open(args.out2, 'wb', compresslevel=args.compression_level)


    manager = Manager()
    read_counter = manager.Value('i', 0)
    reads_removed = manager.Value('i', 0)
    reads_encountered = manager.dict()
    
    rm_dup = partial(remove_duplicates, read_counter, reads_removed, reads_encountered)

    with manager.Pool(3) as wp:
        for read_1, read_2 in wp.imap(rm_dup, (iter(fq1), iter(fq2))):
            out_1.write((read_1 + '\n').encode())
            out_2.write((read_2 + '\n').encode())
    
    out_1.close()
    out_2.close()

    
    

          
    
if __name__ == '__main__':
    main()

