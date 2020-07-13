#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith

Script generates a bed file of restriction fragment locations in a given genome.

"""

import argparse
import sys
import os
import pysam
import re
#import pdb

def get_parser(parser=None):

    if not parser:
        parser= argparse.ArgumentParser()

    parser.add_argument('-i', '--input_fasta', help='fasta file to parse', required=True)
    parser.add_argument('-o', '--output_file', help='output file name',
                   default='digested.bed')

    enzyme_group = parser.add_mutually_exclusive_group(required=True)
    enzyme_group.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
                              default='DpnII')
    enzyme_group.add_argument('-s', '--cut_sequence', help='Sequence of restriction site',
                              default='GATC')
    parser.add_argument('-l', '--logfile', help='filename for logfile',
                        default='test.log')
    return parser


def get_re_site(cut_sequence=None,
                restriction_enzyme=None):
    known_enzymes = {'dpnii': 'GATC',
                     'mboi': 'GATC',
                     'hindiii': 'AAGCTT',
                     'ecori': 'GAATTC'}
    if cut_sequence:
        return cut_sequence.upper()
    elif restriction_enzyme.lower() in known_enzymes:
        return known_enzymes.get(restriction_enzyme.lower())
    else:
        raise ValueError('No restriction site or recognised enzyme provided')

def main(input_fasta,
         restriction_enzyme=None,
         cut_sequence=None,
         logfile=None,
         output_file=None,
         cut_offset=0):

    cut_sequence = get_re_site(restriction_enzyme=restriction_enzyme,
                               cut_sequence=cut_sequence)

    with open(logfile, 'w') as log,\
         open(output_file, 'w') as bed_out,\
         pysam.FastxFile(input_fasta) as fasta_file:

        for seq_entry in fasta_file:

            # Find match positions of the restriction enzyme sequence
            seq_length = len(seq_entry.sequence)
            match_positions = [m.start() for m in re.finditer(cut_sequence, seq_entry.sequence.upper())]

            # iterate through matches and write to bed file
            slice_start = 0
            for match_index, match_pos in enumerate(match_positions):
                slice_end = match_pos + cut_offset

                if slice_start != slice_end:
                    slice_name = f'{restriction_enzyme}_{seq_entry.name}_{match_index}'
                    bed_out.write(f'{seq_entry.name}\t{slice_start}\t{slice_end}\t{slice_name}\n')
                slice_start = slice_end

            # handle last slice
            if slice_start != seq_length:
                slice_end = seq_length
                slice_name = f'{restriction_enzyme}_{seq_entry.name}_{match_index}'
                bed_out.write(f'{seq_entry.name}\t{slice_start}\t{slice_end}\t{slice_name}\n')

            # Print total slices per chr to log file
            log.write(f'{seq_entry.name}: {len(match_positions)}\n')


if __name__ == '__main__':
    main(**vars(get_parser().parse_args()))
