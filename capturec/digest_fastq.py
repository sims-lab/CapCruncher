#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 19:36:52 2020

@author: asmith
"""
import argparse
import multiprocessing as mp
import os
import re
import sys
from collections import Counter
import time

import pandas as pd
import pysam
from pysam import FastxFile
from xopen import xopen

# Make sure the script can find capturec scripts
SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(PACKAGE_DIR)

from .digest_genome import get_re_site

def parse_args():
    parser = argparse.ArgumentParser(prog='digest_fastq')

    subparsers = parser.add_subparsers(
        help='Run in either flashed or unflashed',
        dest='command')
    parser_flashed = subparsers.add_parser('flashed', help='For flashed reads')
    parser_flashed.add_argument('-i', '--input_fastq',
                                help='fastq file to parse',
                                required=True)

    parser_unflashed = subparsers.add_parser(
        'unflashed', help='For unflashed reads')
    parser_unflashed.add_argument('-1', '--fq1',
                                  help='fastq file containing read 1 to parse',
                                  required=True)
    parser_unflashed.add_argument('-2', '--fq2',
                                  help='fastq file containing read 2 to parse',
                                  required=True)

    for subparser in [parser_unflashed, parser_flashed]:
        subparser = add_shared_options(subparser)

    return parser.parse_args()


def add_shared_options(parser):
    parser.add_argument('-o', '--output_file', help='output file name',
                        default='digested.fastq.gz')

    enzyme_group = parser.add_mutually_exclusive_group(required=True)
    enzyme_group.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
                              default='DpnII')
    enzyme_group.add_argument('-s', '--cut_sequence',
                               help='Sequence of restriction site')

    parser.add_argument('-m', '--minimum_slice_length', help='Shortest length for a slice to be output',
                        default=20, type=int)
    parser.add_argument('--stats_file',
                        help='stats_file_prefix', default='stats.log')
    parser.add_argument('-c', '--compression_level', help='Level of gzip compression (1-9 with 9 being the most compressed/slowest)',
                        default=6, type=int)
    parser.add_argument('--keep_cutsite', help='Determines if cutsite is stripped from the start of each slice',
                        action='store_true', default=False)
    parser.add_argument(
        '--buffer', help='Number of reads to process before writing output', default=10000, type=int)

    parser.add_argument(
        '-p', '--n_digestion_processes',help='Number of digestion processes to spawn', default=1, type=int)

    return parser

class DigestedRead():
    '''Class performs in silico digestion of fastq reads and contains relevant stats.

        Args:
            read - read in fastq format to be digested.
            cutsite - compiled regex for the site of digestion.
            flashed - determines if reads have been combined (i.e. using FLASh).
            minimum_slice_length - determines the minimum number of basepairs
                                   for a valid slice.
            slice_offset - all output slices will be labled by slice_offset
                           and slice number. (useful for unflashed pairs)
            keep_cutsite - removes the cutsite from the output slice if False

        Attributes:
            recognition_sites - list of identified restriction sites in the read
            slices_total_counter - number of un-validated slices in the read
            slices_valid_counter - number of validated slices in the read
            slices - list of slices in fastq format (string)

        '''

    def __init__(self,
                 read: pysam.FastqProxy,
                 cutsite: re.compile,
                 flashed=False,
                 minimum_slice_length=0,
                 slice_offset=0,
                 keep_cutsite=False):

        self.read = read  # object with attributes: name, sequence, quality
        self.cutsite = cutsite  # compiled regex designating the cutsite
        self.min_slice_len = minimum_slice_length
        self.flashed = flashed
        self.read_type = 'flashed' if flashed else 'unflashed'
        # Determines if the recognition site is removed from each slice
        self.keep_cutsite = keep_cutsite

        self.recognition_sites = ([site.start() for site in cutsite.finditer(self.read.sequence.upper())] +
                                  [len(self.read.sequence)])  # Find the start location of each recognition site
                                                              # the end position of the sequence is also added
        self.slices_total_counter = 0
        self.slices_valid_counter = 0
        self.slice_offset = slice_offset  # Enables adjusting the slice output number

        self.slices = self.get_slices()

    def get_slices(self):
        '''Splits the read sequence into slices based on the presence of a
           restriction site.
        '''
        slice_start = 0
        slices_lst = []
        for site in self.recognition_sites:

            self.slices_total_counter += 1
            slice_end = site
            slice_length = slice_end - slice_start

            if not self.keep_cutsite:  # Check if the cutsite needs to be removed
                cutsite_removed = re.sub(
                    f'^{self.cutsite.pattern}', '', self.read.sequence[slice_start:slice_end])
                slice_shift = len(
                    self.read.sequence[slice_start:slice_end]) - len(cutsite_removed)
                slice_start += slice_shift  # Shift the slice by the length of the removed cutsite
                slice_length = slice_end - slice_start

            # Make a temporary variable to hold the unvalidated slice
            s = '\n'.join([f'@{self.read.name}|{self.read_type}|{self.slices_valid_counter + self.slice_offset}',
                           self.read.sequence[slice_start:slice_end],
                           '+',
                           self.read.quality[slice_start:slice_end]])

            if slice_length >= self.min_slice_len:  # Confirm that the slice meets minimum length requirement

                # Only allow a slice to be recorded if the read is unflashed or digestion has occured
                if (not self.flashed) or (self.flashed and slice_length < len(self.read.sequence)):
                    self.slices_valid_counter += 1
                    slices_lst.append(s)

            slice_start = slice_end

        return slices_lst

    def __str__(self):
        if self.slices:
            return '\n'.join(self.slices) + '\n'
        else:
            return ''

def read_fastq(fq,
               outq,
               n_workers=1,
               buffer=10000):
    '''Reads fastq file and places reads into a queue'''
    r_buffer = []
    for rc, read in enumerate(fq):
        r_buffer.append(read)

        if rc % buffer == 0 and not rc == 0:
            outq.put(r_buffer)
            r_buffer = []
            print(f'Processed {rc} reads')

    outq.put(r_buffer) # Add the reads that don't fit the buffer size
    for _ in range(n_workers):
        outq.put('TER') # Places a terminator in the queue for every spawned worker


def read_paired_fastqs(fq1,
                       fq2,
                       outq,
                       n_workers=1,
                       buffer=10000):
    '''Reads R1 and R2 fastq files and places paired reads into a queue'''
    r_buffer = []
    for rc, (r1, r2) in enumerate(zip(fq1, fq2)):
        r_buffer.append((r1, r2))

        if rc % buffer == 0 and not rc == 0:
            outq.put(r_buffer)
            r_buffer = []
            print(f'Processed {rc} reads')

    outq.put(r_buffer)
    for _ in range(n_workers):
        outq.put('TER') # Places a terminator in the queue for every spawned worker



def digest_read_flashed(inq, outq, statq, **kwargs):

    read_buffer = []
    stat_buffer = []
    reads = inq.get()

    while not reads == 'TER':
        for read in reads:
            sliced_read = DigestedRead(read, **kwargs)
            sliced_read_str = str(sliced_read)
            if sliced_read_str: # Only append if valid slices present
                read_buffer.append(sliced_read_str)

            stat_buffer.append((sliced_read.slices_total_counter,
                                sliced_read.slices_valid_counter,
                                )
                              )

        outq.put(read_buffer) # Add list of digested reads to the queue
        statq.put(stat_buffer)
        read_buffer = [] # Clear buffer
        stat_buffer = [] # Clear buffer
        reads = inq.get() # Get new list of reads

    outq.put('TER')
    statq.put('TER')



def digest_read_unflashed(inq, outq, statq, **kwargs):

    read_buffer = []
    stat_buffer = []
    reads = inq.get()
    while not reads == 'TER': # Checks to see if the queue has been terminated
        for r1, r2 in reads:
            sliced_read_1 = DigestedRead(r1, **kwargs) # Digest read 1
            kwargs['slice_offset'] = sliced_read_1.slices_valid_counter # Update slice offset
            sliced_read_2 = DigestedRead(r2, **kwargs) # Digest read 2
            kwargs['slice_offset'] = 0 # Reset slice offset

            s1, s2 = str(sliced_read_1), str(sliced_read_2)

            if s1 and s2: # Only store of both reads have valid slices
                read_buffer.append(f'{s1}{s2}')

            stat_buffer.append((sliced_read_1.slices_total_counter,
                                sliced_read_1.slices_valid_counter,
                                sliced_read_2.slices_total_counter,
                                sliced_read_2.slices_valid_counter,
                                    ))



        outq.put(read_buffer)
        statq.put(stat_buffer)
        read_buffer = []
        stat_buffer = []
        reads = inq.get()

    outq.put('TER')
    statq.put('TER')


def write_to_fastq(inq,
                   outfile='out.fq.gz',
                   compression_level=5):

    '''Writes all digested reads to the appropriate file'''
    with xopen(filename=outfile, mode='wb', compresslevel=compression_level) as f:

        reads = inq.get()
        while not reads == 'TER':
            f.write(''.join(reads).encode())
            reads = inq.get()

def collate_stats_flashed(inq, stats_file='out.log'):

    counts = inq.get()
    total_counter = Counter()
    valid_counter = Counter()
    while not counts == 'TER':
        total, valid = zip(*counts)
        total_counter = total_counter +  Counter(total)
        valid_counter = valid_counter + Counter(valid)
        counts = inq.get()


    dframes = []
    for key, values in zip(['total', 'valid'],
                           [total_counter, valid_counter]):

        dframes.append(pd.DataFrame.from_dict(values, orient='index')
                        .reset_index()
                        .rename(columns={'index':'bin', 0: 'frequency'})
                        .assign(stat=key))

    df_stats = pd.concat(dframes)
    df_stats['read_type'] = 'flashed'
    df_stats['bin'] = df_stats['bin'].astype(int)

    (df_stats.sort_values('bin')
             .to_csv(stats_file, index=False))

def collate_stats_unflashed(inq, stats_file='out.log'):

    counts = inq.get()
    counters = {'r1': {'total': Counter(),
                       'valid': Counter()
                       },
                'r2': {'total': Counter(),
                        'valid': Counter()}
                }
    while not counts == 'TER':
        r1_total, r1_valid, r2_total, r2_valid = zip(*counts)
        counters['r1']['total'] = counters['r1']['total'] + Counter(r1_total)
        counters['r1']['valid'] = counters['r1']['valid'] + Counter(r1_valid)
        counters['r2']['total'] = counters['r2']['total'] + Counter(r2_total)
        counters['r2']['valid'] = counters['r2']['valid'] + Counter(r2_valid)
        counts = inq.get()

    dframes = []
    for read in counters:
        for key, values in counters[read].items():

            dframes.append(pd.DataFrame.from_dict(values, orient='index')
                            .reset_index()
                            .rename(columns={'index':'bin', 0: 'frequency'})
                            .assign(stat=key, read_type=read))

    df_stats = pd.concat(dframes)
    df_stats['bin'] = df_stats['bin'].astype(int)

    (df_stats.sort_values('bin')
             .to_csv(stats_file, index=False))


def main(command,
         input_fastq=None,
         fq1=None,
         fq2=None,
         cut_sequence=None,
         restriction_enzyme=None,
         keep_cutsite=False,
         output_file='out.fastq.gz',
         minimum_slice_length=18,
         stats_file=None,
         compression_level=5,
         n_digestion_processes=1,
         buffer=10000):


    # Set up multiprocessing variables
    inputq = mp.SimpleQueue() # reads are placed into this queue for processing
    writeq = mp.SimpleQueue() # digested reads are placed into the queue for writing
    statq = mp.SimpleQueue() # stats queue

    # Variables
    cut_site = re.compile(get_re_site(cut_sequence, restriction_enzyme))

    # # Remove stats_file if present
    # if not stats_file == sys.stdout and os.path.exists(stats_file):
    #     os.remove(stats_file)


    if command == 'flashed':  # Checks the subcommand to see in which mode to run

        fq = FastxFile(input_fastq)

        # Writer process
        writer = mp.Process(target=write_to_fastq, args=(writeq,))
        writer.start()

        # Define processes
        processes = [mp.Process(target=read_fastq,
                                args=(fq, inputq),
                                kwargs={'n_workers': n_digestion_processes}
                                ),
                     mp.Process(target=collate_stats_flashed,
                                args=(statq,),
                                kwargs={'stats_file': stats_file},
                                )
                    ]

        processes_repeated = [mp.Process(target=digest_read_flashed,
                                         args=(inputq, writeq, statq),
                                         kwargs={'cutsite': cut_site,
                                                 'flashed': True,
                                                 'minimum_slice_length':
                                                  minimum_slice_length,
                                                 'keep_cutsite': keep_cutsite})
                              for i in range(n_digestion_processes)]

        processes = processes + processes_repeated

    elif command == 'unflashed':

        fq1, fq2 = FastxFile(fq1), FastxFile(fq2)

        # Writer process
        writer = mp.Process(target=write_to_fastq, args=(writeq,))
        writer.start()

        # Define processes
        processes = [mp.Process(target=read_paired_fastqs,
                                args=(fq1, fq2, inputq),
                                kwargs={'n_workers': n_digestion_processes}
                                ),
                     mp.Process(target=collate_stats_unflashed,
                                args=(statq,),
                                kwargs={'stats_file': stats_file},
                                )
                    ]

        processes_repeated = [mp.Process(target=digest_read_unflashed,
                                         args=(inputq, writeq, statq),
                                         kwargs={'cutsite': cut_site,
                                                 'flashed': False,
                                                 'minimum_slice_length': min_slice_len,
                                                 'keep_cutsite': keep_cutsite,
                                                 'slice_offset': 0})
                              for i in range(n_digestion_processes)]

        processes = processes + processes_repeated

    # Start all processes
    for proc in processes:
        proc.start()

    # Join processes (wait for processes to finish before the main process)
    writer.join()

    for proc in processes:
        proc.join(10)
        proc.terminate()


if __name__ == '__main__':
    main(**vars(parse_args()))
