#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""

import argparse
import sys
import os
import pysam
import re
#import pdb

p = argparse.ArgumentParser()
p.add_argument('-i', '--input_fasta', help='fasta file to parse')
p.add_argument('-o', '--output_file', help='output file name', 
               default='digested.bed')
p.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
               default='DpnII')
p.add_argument('-s', '--cut_sequence', help='Sequence of restriction site',
               default='GATC')
p.add_argument('-l', '--logfile', help='filename for logfile',
               default=sys.stdout)
args = p.parse_args()

cut_offset = 0

# command line input: -i test.fa -l fasta.log

# assertions - check input file exists
assert os.path.isfile(args.input_fasta), "Input fasta file not found"

def main():
    logfile = open(args.logfile, 'w')
    with open(args.output_file, 'w') as bed_out:
        with pysam.FastxFile(args.input_fasta) as fasta_file:
            for seq_entry in fasta_file:
                # Find match positions of the restriction enzyme sequence
                seq_length = len(seq_entry.sequence)
                match_positions = [m.start() for m in re.finditer(args.cut_sequence, seq_entry.sequence.upper())]
                # iterate through matches and write to bed file
                slice_start = 0
                for match_index, match_pos in enumerate(match_positions):
                    slice_end = match_pos + cut_offset
                    if slice_start != slice_end:
                        slice_name = f'{args.restriction_enzyme}_{seq_entry.name}_{match_index}'
                        bed_out.write(f'{seq_entry.name}\t{slice_start}\t{slice_end}\t{slice_name}\n')
                    slice_start = slice_end
                # handle last slice
                if slice_start != seq_length:
                    slice_end = seq_length
                    slice_name = f'{args.restriction_enzyme}_{seq_entry.name}_{match_index}'
                    bed_out.write(f'{seq_entry.name}\t{slice_start}\t{slice_end}\t{slice_name}\n')
                # Print total slices per chr to log file
                logfile.write(f'{seq_entry.name}: {len(match_positions)}\n')
    logfile.close()

if __name__ == '__main__':
    main()