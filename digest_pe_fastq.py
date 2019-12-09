import argparse
import itertools
import os
import re
import sys
import gzip
from collections import defaultdict
import numba
from multiprocessing import Pool
import time

from pysam import FastxFile

p = argparse.ArgumentParser()
p.add_argument('-1', '--fq1', help='fastq file containing read 1 to parse')
p.add_argument('-2', '--fq2', help='fastq file containing read 2 to parse')
p.add_argument('-o', '--output_file', help='output file name', 
               default='digested_PE.fastq.gz')
p.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
               default='DpnII')
p.add_argument('-s', '--cut_sequence', help='Sequence of restriction site', default=None)
p.add_argument('-m', '--minimum_slice_length', help='Shortest length for a slice to be output',
               default=0, type=int)
p.add_argument('-l', '--logfile', help='filename for logfile',
               default=sys.stdout)
p.add_argument('-c', '--compression_level', help='Level of gzip compression (1-9 with 9 being the most compressed/slowest)',
               default=6, type=int)

args = p.parse_args()


# assertions - check input file exists
assert os.path.isfile(args.fq1), "Input fastq file containing read 1 not found"
assert os.path.isfile(args.fq2), "Input fastq file containing read 2 not found"


known_enzymes = {'dpnii': 'GATC',
                 'mboi': 'GATC',}

if args.cut_sequence != None:
    cut_sequence = args.cut_sequence
elif args.restriction_enzyme:
    cut_sequence = known_enzymes.get(args.restriction_enzyme.lower())


def find_cut_sites(seq):
    return [m.start() for m in re.finditer(cut_sequence, seq)] + [len(seq)]

def get_slices(seq, match_positions, min_slice_len):
    
    if len(match_positions) > 0:
        slice_start = 0
        for match_pos in match_positions:
            slice_data = dict()
            slice_end = match_pos
            slice_length = slice_end - slice_start
            
            # Make sure that the slice meets the minimum length requirement
            if slice_length > min_slice_len and slice_end != 0:
                slice_data['sequence'] = seq[slice_start:slice_end]
                slice_data['slice_start'] = slice_start
                slice_data['slice_end'] = slice_end            
                yield slice_data
                slice_start = slice_end
            else:
            # Return None if the slice does not pass the threshold
                yield None
                slice_start = slice_end 

    # In this case the cutsite might not be observed so will keep any undigested reads
    else:
        slice_data = dict()
        slice_data['sequence'] = seq
        slice_data['slice_start'] = 0
        slice_data['slice_end'] = len(seq)
        yield slice_data


def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return open(fn, 'w')
    else:
        return fn

def main():

    r1_cut_site_counts = defaultdict(int)
    r2_cut_site_counts = defaultdict(int)
    min_slice_len = args.minimum_slice_length
    invalid_slice_counter = 0

    with gzip.open(args.output_file, 'wb', compresslevel=args.compression_level) as fastq_out:
        
        fq1, fq2 = FastxFile(args.fq1), FastxFile(args.fq2)
        
        #Iterate over both fq files
        for seq_counter, (r1, r2) in enumerate(zip(fq1,fq2)):

            if seq_counter % 10000 == 0:
                print(f'Processed {seq_counter} reads')

            #Find cut sites
            r1_match_positions, r2_match_positions = [find_cut_sites(seq) for seq in  [r1.sequence, r2.sequence]]

            #Record the number of cuts observed in read 1 and read 2 (-1 offset to adjust for recording the full read)
            r1_cut_site_counts[len(r1_match_positions) - 1] += 1
            r2_cut_site_counts[len(r2_match_positions) - 1] += 1


            #If insilico digestion has occured then split the read into slices
            r1_slices, r2_slices = [get_slices(r_sequence, r_match_positions, min_slice_len) 
                                    for r_sequence, r_match_positions in
                                    [(r1.sequence, r1_match_positions), (r2.sequence, r2_match_positions),]
                                    ]
  
              # Write the slices for read 1 to a new fastq file
            for ii, s in enumerate(r1_slices):
                if s:
                    fastq_out.write(f'@{r1.name}|PE1|{ii}\n'.encode())
                    fastq_out.write(f'{s["sequence"]}\n'.encode())
                    fastq_out.write('+\n'.encode())
                    fastq_out.write(f'{r1.quality[s["slice_start"]:s["slice_end"]]}\n'.encode())
                else:
                    # If a slice does not pass the threshold then count these
                    invalid_slice_counter += 1

            # Repeat for read 2
            for ii, s in enumerate(r2_slices):
                if s:
                    fastq_out.write(f'@{r2.name}|PE2|{ii}\n'.encode())
                    fastq_out.write(f'{s["sequence"]}\n'.encode())
                    fastq_out.write('+\n'.encode())
                    fastq_out.write(f'{r2.quality[s["slice_start"]:s["slice_end"]]}\n'.encode())
                else:
                    invalid_slice_counter += 1
            

    with open_logfile(args.logfile) as logfile:
        print('==============')
        logfile.write(f'Number of read pairs processed {seq_counter}\n')
        
        logfile.write('Frequency of slices observed in read 1:\n')
        for k in sorted(r1_cut_site_counts):
            logfile.write(f'{k}: {r1_cut_site_counts[k]}\n')
        
        logfile.write('Frequency of slices observed in read 2:\n')
        for k in sorted(r2_cut_site_counts):
            logfile.write(f'{k}: {r2_cut_site_counts[k]}\n')
          
        logfile.write(f'Number of slices below the minimum slice length: {invalid_slice_counter}\n')
        print('==============')

if __name__ == '__main__':
    main()