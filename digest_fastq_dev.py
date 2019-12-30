import argparse
import gzip
import itertools
import os
import re
import sys
import time
from collections import defaultdict

from pysam import FastxFile

parser = argparse.ArgumentParser(prog='digest_fastq')

subparsers = parser.add_subparsers(help='Run in either flashed or unflashed', dest='command')
parser_flashed = subparsers.add_parser('flashed', help='For flashed reads')
parser_flashed.add_argument('-i', '--input_fastq', help='fastq file to parse')

parser_unflashed = subparsers.add_parser('unflashed', help='For unflashed reads')
parser_unflashed.add_argument('-1', '--fq1', help='fastq file containing read 1 to parse')
parser_unflashed.add_argument('-2', '--fq2', help='fastq file containing read 2 to parse')

parser.add_argument('-o', '--output_file', help='output file name', 
               default='digested.fastq.gz')
parser.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
               default='DpnII')
parser.add_argument('-s', '--cut_sequence', help='Sequence of restriction site')
parser.add_argument('-m', '--minimum_slice_length', help='Shortest length for a slice to be output',
               default=0, type=int)
parser.add_argument('-l', '--logfile', help='filename for logfile',
               default=sys.stdout)
parser.add_argument('-c', '--compression_level', help='Level of gzip compression (1-9 with 9 being the most compressed/slowest)',
               default=6, type=int)

args = parser.parse_args()


class DigestedRead():
    def __init__(self, read, cutsite, flashed=False, minimum_slice_length=0, slice_offset=0):
        self.read = read
        self.cutsite = cutsite
        self.min_slice_len = minimum_slice_length
        self.flashed = flashed
        self.read_type = 'flashed' if flashed else 'unflashed'

        self.recognition_sites = [site.start() for site in cutsite.finditer(self.read.sequence)] + [len(self.read.sequence)]
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

            if slice_length > self.min_slice_len and slice_end != 0:
                s = '\n'.join([f'{self.read.name}|{self.read_type}|{self.slices_valid_counter + self.slice_offset}',
                               self.read.sequence[slice_start:slice_end],
                               '+',
                               self.read.quality[slice_start:slice_end]])
                
                if (not self.flashed) or (self.flashed and slice_length != len(self.read.sequence)):
                    self.slices_valid_counter += 1
                    slices_lst.append(s)

        return slices_lst
    
    def get_slices_string(self):
        if self.slices:
            return '\n'.join(self.slices)


def open_logfile(fn):
    if fn == sys.stdout:
        return fn
    else:
        return open(fn, 'w')

def get_re_site(cut_sequence=None, restriction_enzyme=None):
    known_enzymes = {'dpnii': 'GATC',
                     'mboi': 'GATC',}
    
    if cut_sequence:
        return cut_sequence
    elif restriction_enzyme:
        return known_enzymes.get(args.restriction_enzyme.lower())
    else:
        raise ValueError('No restriction site or recognised enzyme provided')


def main():
    
    slice_total_dict = {'flashed': defaultdict(int),
                        'read_1': defaultdict(int),
                        'read_2': defaultdict(int),
                        }
    
    slice_valid_dict = {'flashed': defaultdict(int),
                        'read_1': defaultdict(int),
                        'read_2': defaultdict(int),
                        }
    
    min_slice_len = args.minimum_slice_length
    cut_site = re.compile(get_re_site(args.cut_sequence, args.restriction_enzyme))
    

    with gzip.open(args.output_file, 'w', compresslevel=args.compression_level) as fastq_out: 
        
        if args.command == 'flashed':    
            for seq_counter, read in enumerate(FastxFile(args.input_fastq)):

                if seq_counter % 10000 == 0:
                    print(f'Processed {seq_counter} reads')

                sliced_read = DigestedRead(read,
                                           cut_site,
                                           flashed=True,
                                           minimum_slice_length=min_slice_len)
                
                if sliced_read.slices_string:
                    s = (sliced_read.slices_string + '\n').encode()
                    fastq_out.write(s)
                
                slice_total_dict['flashed'][sliced_read.slices_total_counter] += 1
                slice_valid_dict['flashed'][sliced_read.slices_valid_counter] += 1
        
        elif args.command == 'unflashed':
            for seq_counter, (read1, read2) in enumerate(zip(FastxFile(args.fq1), FastxFile(args.fq2))):
                
                if seq_counter % 10000 == 0:
                    print(f'Processed {seq_counter} reads')
                

                sliced_read_1 = DigestedRead(read1,
                                             cut_site,
                                             flashed=False,
                                             minimum_slice_length=min_slice_len)

               
                sliced_read_2 = DigestedRead(read2,
                                             cut_site,
                                             flashed=False,
                                             minimum_slice_length=min_slice_len,
                                             slice_offset=sliced_read_1.slices_valid_counter)
                

                for read_name, sliced_read in zip(['read_1', 'read_2'], [sliced_read_1, sliced_read_2]):
                    fastq_out.write((sliced_read.slices_string + '\n').encode())
                    slice_total_dict[read_name][sliced_read.slices_total_counter] += 1
                    slice_valid_dict[read_name][sliced_read.slices_valid_counter] += 1
        

      
        logfile = open_logfile(args.logfile)
        
        logfile.write(f'Reads processed\t{seq_counter}\n')
        for read_type in slice_total_dict:
            total_slices = sum(k * v for k, v in slice_total_dict[read_type].items())
            valid_slices = sum(k * v for k, v in slice_valid_dict[read_type].items())

            if total_slices:
                logfile.write(f'Total slices {read_type}\t{total_slices}\n')
                logfile.write(f'Valid slices {read_type}\t{valid_slices}\n')
                logfile.write(f'Valid slices for {read_type} per read:\n')
                
                for counter in sorted(slice_valid_dict[read_type]):                  
                    logfile.write(f'{counter}\t{slice_valid_dict[read_type][counter]}\n')

        logfile.close()


  




if __name__ == "__main__":
    main()