#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:45:09 2020

@author: asmith

Script splits a fastq into specified chunks
"""
import argparse
import os
import sys
from xopen import xopen
from pysam import FastxFile

def get_parser(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_fastq', help='BAM file to parse')
    parser.add_argument('--chunksize', help='Number of reads per output file', default=1000000, type=int)
    parser.add_argument('-c', '--compression_level', help='Level of gzip compression (1-9 with 9 being the most compressed/slowest)',
                        default=6, type=int)
    parser.add_argument('-n','--output_prefix', help='output prefix')

    return parser

def main(input_fastq,
         output_prefix,
         compression_level=5,
         chunksize=1000000):

    split_counter = 0
    for read_counter, read in enumerate(FastxFile(input_fastq)):

        processed_read = f'{read}\n'.encode()

        if read_counter == 0:
            out_name = f'{output_prefix}_{split_counter}.fastq.gz'
            out_handle = xopen(out_name,
                               'wb',
                               compresslevel=compression_level,
                               threads=8)
            out_handle.write(processed_read)

        elif read_counter % chunksize == 0:
            split_counter += 1
            out_handle.close()
            out_name = f'{output_prefix}_{split_counter}.fastq.gz'
            out_handle =  xopen(out_name,
                               'wb',
                               compresslevel=compression_level,
                               threads=8)
            out_handle.write(processed_read)
        else:
            out_handle.write(processed_read)

    out_handle.close()

if __name__ == '__main__':
    args = get_parser().parse_args()
    main(**vars(args))
