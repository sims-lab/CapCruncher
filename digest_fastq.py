#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith, dsims
"""

import argparse
import sys
import os
import pysam
import gzip
import re
#import pdb

p = argparse.ArgumentParser()
p.add_argument('-i', '--input_fastq', help='fastq file to parse')
p.add_argument('-o', '--output_file', help='output file name', 
               default='digested.fastq.gz')
p.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
               default='DpnII')
p.add_argument('-s', '--cut_sequence', help='Sequence of restriction site')
p.add_argument('-m', '--minimum_slice_length', help='Shortest length for a slice to be output',
               default='0')
p.add_argument('-l', '--logfile', help='filename for logfile',
               default=sys.stdout)
args = p.parse_args()

# assertions - check input file exists
assert os.path.isfile(args.input_fastq), "Input fastq file not found"

# Set cut sequence for known cutters
cut_sequence = None
if args.cut_sequence != None:
    cut_sequence = args.cut_sequence
elif args.restriction_enzyme.lower() == 'dpnii':
    cut_sequence = 'GATC'
else:
    # raise error
    pass
min_slice_length = int(args.minimum_slice_length)

# command line options: -i test.fq -l test.log

cut_offset = 0

def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return open(fn, 'w')
    else:
        return fn

def main():
    with gzip.open(args.output_file, 'wb') as fastq_out:
        # Compile regular expression 
        cutsite_counter = 0
        total_slices = 0
        short_slice_counter = 0
        cut_site_counts = dict()
        for seq_counter, seq_entry in enumerate(pysam.FastqFile(args.input_fastq)):
            # Split the sequence using the restriction enzyme sequence
            seq_length = len(seq_entry.sequence)
            match_positions = [m.start() for m in re.finditer(cut_sequence, seq_entry.sequence)]
            
            # Record how many re digestion sites were found in the read
                        
            if len(match_positions) not in cut_site_counts:
                cut_site_counts[len(match_positions)] = 0  
            cut_site_counts[len(match_positions)] += 1

            if len(match_positions) > 0:
                cutsite_counter += 1
                slice_start = 0
                match_index = 0
                for match_pos in match_positions:
                    #current_slice = last_slice + slice_length
                    slice_end = match_pos + cut_offset
                    slice_length = slice_end - slice_start
                    if slice_length > min_slice_length and slice_end != 0:
                        total_slices += 1
                        fastq_out.write(f'@{seq_entry.name}|PE1|{match_index}\n'.encode())
                        fastq_out.write(f'{seq_entry.sequence[slice_start:slice_end]}\n'.encode())
                        fastq_out.write('+\n'.encode())
                        fastq_out.write(f'{seq_entry.quality[slice_start:slice_end]}\n'.encode())
                        match_index +=1
                    else:
                        short_slice_counter += 1
                    slice_start = slice_end
                # handle last slice
                if slice_start != seq_length:
                    match_index += 1
                    slice_end = seq_length
                    slice_length = slice_end - slice_start
                    if slice_length > min_slice_length:
                        total_slices += 1
                        fastq_out.write(f'@{seq_entry.name}|PE1|{match_index}\n'.encode())
                        fastq_out.write(f'{seq_entry.sequence[slice_start:slice_end]}\n'.encode())
                        fastq_out.write('+\n'.encode())
                        fastq_out.write(f'{seq_entry.quality[slice_start:slice_end]}\n'.encode())

    with open_logfile(args.logfile) as logfile:
        logfile.write(f'Records processed: {seq_counter+1}\n')
        logfile.write(f'Records with cutsites: {cutsite_counter}\n')
        logfile.write(f'slices output: {total_slices}\n')
        logfile.write(f'Too short slices: {short_slice_counter}\n')
        logfile.write(f'Cut count histogram:\n')
        for cut_no in sorted(cut_site_counts.keys()):
            logfile.write(f'{cut_no}\t{cut_site_counts[cut_no]}\n')


if __name__ == '__main__':
    main()