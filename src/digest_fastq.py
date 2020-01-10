import argparse
import gzip
import os
import re
import sys
import pandas as pd
from collections import defaultdict

from pysam import FastxFile

parser = argparse.ArgumentParser(prog='digest_fastq')
parser.add_argument('-o', '--output_file', help='output file name',
                    default='digested.fastq.gz')
parser.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
                    default='DpnII')
parser.add_argument('-s', '--cut_sequence',
                    help='Sequence of restriction site')
parser.add_argument('-m', '--minimum_slice_length', help='Shortest length for a slice to be output',
                    default=0, type=int)
parser.add_argument('-l', '--logfile',
                    help='logfile_prefix', default=sys.stdout)
parser.add_argument('-c', '--compression_level', help='Level of gzip compression (1-9 with 9 being the most compressed/slowest)',
                    default=6, type=int)
parser.add_argument('--keep_cutsite', help='Determines if cutsite is stripped from the start of each slice',
                    action='store_true', default=False)

subparsers = parser.add_subparsers(
    help='Run in either flashed or unflashed', dest='command')
parser_flashed = subparsers.add_parser('flashed', help='For flashed reads')
parser_flashed.add_argument('-i', '--input_fastq', help='fastq file to parse')

parser_unflashed = subparsers.add_parser(
    'unflashed', help='For unflashed reads')
parser_unflashed.add_argument(
    '-1', '--fq1', help='fastq file containing read 1 to parse')
parser_unflashed.add_argument(
    '-2', '--fq2', help='fastq file containing read 2 to parse')


args = parser.parse_args()


class DigestedRead():
    def __init__(self, read, cutsite, flashed=False, minimum_slice_length=0, slice_offset=0, keep_cutsite=False):
        self.read = read
        self.cutsite = cutsite
        self.min_slice_len = minimum_slice_length
        self.flashed = flashed
        self.read_type = 'flashed' if flashed else 'unflashed'
        self.keep_cutsite = keep_cutsite

        self.recognition_sites = ([site.start() for site in cutsite.finditer(self.read.sequence.upper())] + 
                                  [len(self.read.sequence)])
        self.slices_total_counter = 0
        self.slices_valid_counter = 0
        self.slice_offset = slice_offset

        self.slices = self.get_slices()
        self.slices_string = self.get_slices_string()

    def get_slices(self):
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

            s = '\n'.join([f'@{self.read.name}|{self.read_type}|{self.slices_valid_counter + self.slice_offset}',
                           self.read.sequence[slice_start:slice_end],
                           '+',
                           self.read.quality[slice_start:slice_end]])

            if slice_length >= self.min_slice_len:
                if (not self.flashed) or (self.flashed and slice_length < len(self.read.sequence)):
                    # Only allow a slice to be recorded if the read is unflashed or digestion has occured
                    self.slices_valid_counter += 1
                    slices_lst.append(s)

            slice_start = slice_end

        return slices_lst

    def get_slices_string(self):
        if self.slices:
            return '\n'.join(self.slices)


def get_re_site(cut_sequence=None, restriction_enzyme=None):
    known_enzymes = {'dpnii': 'GATC',
                     'mboi': 'GATC', }

    if cut_sequence:
        return cut_sequence
    elif restriction_enzyme:
        return known_enzymes.get(args.restriction_enzyme.lower())
    else:
        raise ValueError('No restriction site or recognised enzyme provided')


def get_digestion_stats(n_processed, total_slices, valid_slices):
    stats_combined = dict()
    for read_type in total_slices:
        total_count = sum(k * v for k, v in total_slices[read_type].items())

        if total_count:
            stats = {'total_read_pairs_processed': n_processed,
                     'total_slices': sum(k * v for k, v in total_slices[read_type].items()),
                     'total_valid_slices': sum(k * v for k, v in valid_slices[read_type].items())}

            hist = {k: valid_slices[read_type][k] 
                    for k in sorted(valid_slices[read_type])}
            stats.update(hist)
        else:
            stats = {'total_read_pairs_processed': 0,
                        'total_slices': 0,
                        'total_valid_slices': 0}

        stats_combined[read_type] = stats
    
    return pd.DataFrame(stats_combined)

def main():

    total_slices = {'flashed': defaultdict(int),
                        'read_1': defaultdict(int),
                        'read_2': defaultdict(int),
                        }

    valid_slices = {'flashed': defaultdict(int),
                        'read_1': defaultdict(int),
                        'read_2': defaultdict(int),
                        }

    min_slice_len = args.minimum_slice_length
    cut_site = re.compile(get_re_site(args.cut_sequence, args.restriction_enzyme))
    keep_cutsite = args.keep_cutsite

    with gzip.open(args.output_file, 'w', compresslevel=args.compression_level) as fastq_out:

        if args.command == 'flashed':
            for seq_counter, read in enumerate(FastxFile(args.input_fastq)):

                if seq_counter % 10000 == 0:
                    print(f'Processed {seq_counter} reads')

                sliced_read = DigestedRead(read,
                                           cut_site,
                                           flashed=True,
                                           minimum_slice_length=min_slice_len,
                                           keep_cutsite=keep_cutsite)

                if sliced_read.slices_string:
                    s = (sliced_read.slices_string + '\n').encode()
                    fastq_out.write(s)

                total_slices['flashed'][sliced_read.slices_total_counter] += 1
                valid_slices['flashed'][sliced_read.slices_valid_counter] += 1

        elif args.command == 'unflashed':
            for seq_counter, (read1, read2) in enumerate(zip(FastxFile(args.fq1), FastxFile(args.fq2))):

                if seq_counter % 10000 == 0:
                    print(f'Processed {seq_counter} reads')

                sliced_read_1 = DigestedRead(read1,
                                             cut_site,
                                             flashed=False,
                                             minimum_slice_length=min_slice_len,
                                             keep_cutsite=keep_cutsite)

                sliced_read_2 = DigestedRead(read2,
                                             cut_site,
                                             flashed=False,
                                             minimum_slice_length=min_slice_len,
                                             slice_offset=sliced_read_1.slices_valid_counter,
                                             keep_cutsite=keep_cutsite)

                for read_name, sliced_read in zip(['read_1', 'read_2'], [sliced_read_1, sliced_read_2]):
                    if sliced_read.slices_string:
                        s = (sliced_read.slices_string + '\n').encode()
                        fastq_out.write(s)

                    total_slices[read_name][sliced_read.slices_total_counter] += 1
                    valid_slices[read_name][sliced_read.slices_valid_counter] += 1

        
        df_stats = get_digestion_stats(seq_counter, total_slices, valid_slices)
        df_stats.to_csv(f'{args.logfile}.tsv', sep='\t')


if __name__ == "__main__":
    main()
