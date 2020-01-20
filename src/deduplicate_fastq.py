#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""

import argparse
import sys
from pysam import FastxFile
import gzip
import pandas as pd

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
    '''Deals with the logfile being set to stdout'''
    if not isinstance(fn, type(sys.stdout)):
        return open(fn, 'w')
    else:
        return fn

def main():

    fq1, fq2 = FastxFile(args.fq1), FastxFile(args.fq2) # Use pysam to open both fastq files

    with gzip.open(args.out1, 'wb', compresslevel=args.compression_level) as of1,\
         gzip.open(args.out2, 'wb', compresslevel=args.compression_level) as of2:
        
        read_encountered = set() # Used to store all unique read hashes
        
        reads_removed = 0
        for i, (r1, r2) in enumerate(zip(fq1, fq2)):
            read_pair = hash(r1.sequence +  r2.sequence) # Stores the hash of the combined sequence
                                                         # hashing is performed to save memory 
            if i % 100000 == 0:
                print('Processed %s reads' % (i))

            if not read_pair in read_encountered: # If read is unique then write it to a fastq file
                of1.write((str(r1) +'\n').encode())
                of2.write((str(r2) +'\n').encode())
                read_encountered.add(read_pair) # Add the unique read to the set
            else:
                reads_removed += 1
        
        with open_logfile(args.logfile) as logfile:
            logfile.write(f'Read_pairs_processed\t{i+1}\n')
            logfile.write(f'Read_pairs_unique\t{i-reads_removed}\n')
            logfile.write(f'Read_pairs_removed\t{reads_removed}\n')
    
if __name__ == '__main__':
    main()

