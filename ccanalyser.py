#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 21:11:23 2019

@author: davids

TODO


- don't allow multiple reports for the same RE fragment from the same capture (DS)
    - filter step before output
    - integrate fragment into perslice dataframe
- blacklist (DS)
    - ucsc chipseq blacklist as example
- remove duplicate from flashed fastq file or raw fastq (AS)
    - iterate fastq(s) and remove identical reads (both ends identical seq) 
    - use set
    - implement as dedup_fastq.py script
- output per viewpoint - bedgraph
    - groupby pandas
- make trackhub (overlay) (AS)
    - seperate pipeline step
    - colours in yml
- stat file output tab delimited (DS)
    - per fragment
    - per slice
    - PCR duplicate rate - stats
    - Digestion efficiency - capture + exclusion stats
- report / plots (AS)
    - jupyter notebook
    - Jelena Python code?
- unflashed reads  (DS, AS)
    - new digest_pe_fastq.py script (AS)
    - regex
    - read names the same for both ends
    - additional info in reads names - PE1/2, last fragment?
    - iterate both fastq together
    - no digestion stats and digestion site might be missed
    - output single fastq file with all slices together (interleaved)
- tri-c support
- mnase support
- Do we need bam2bed?
- chunk input fastq files and recombining
- wobble is base positions for deduplication
    - do this at merge / final step - will be slow
- more test data

"""

import argparse
import os
import pysam
import pandas as pd
import time
import numba
#import pdb

start = time.time()

# Parse input parameters
p = argparse.ArgumentParser()
p.add_argument('-i', '--input_bam', help='BAM file to parse')
#p.add_argument('-r', '--genome_digest_bed', 
#               help='BED file of restriction enzyme sites in the genome (sorted)')
#p.add_argument('-b', '--blacklist', 
#               help='BED file of blacklist regions')
p.add_argument('-c', '--capture', 
               help='intersection  of reads and capture sites from bedtools')
p.add_argument('-d', '--capture_count',
               help='capture site counts')
p.add_argument('-x', '--exclude', 
               help='intersection  of reads and exclusion sites from bedtools')
p.add_argument('-y', '--exclude_count',
               help='exculsion site counts')
p.add_argument('-o', '--output', help='output file prefix')
p.add_argument('-l', '--logfile', help='log file name')
args = p.parse_args()

chunk_size = 100000
debug = False

# commmand line options: -i test.digest.bam -c miseq.digest.capture.intersect -x miseq.digest.exclude.intersect -d miseq.digest.capture.count.bed -y miseq.digest.exclude.count.bed -o test

# assertions - check all input files exist
assert os.path.isfile(args.input_bam), "Input sam file not found"
#assert os.path.isfile(args.genome_digest_bed), "Genome restriction enzyme digest bed file not found"
#assert os.path.isfile(args.blacklist), "Blacklist bed file not found"
assert os.path.isfile(args.capture), "Capture site overlap file  not found"
assert os.path.isfile(args.capture_count), "Capture site count file  not found"
assert os.path.isfile(args.exclude), "Exclusion site overlap file  not found"
assert os.path.isfile(args.exclude_count), "Exclusion site count file  not found"

# Set default output files
output_prefix = None
if args.output is None:
    output_prefix = os.path.basename(args.input_bam)
else:
    output_prefix = args.output

log_file = None
if args.logfile is None:
    log_file = output_prefix + '.log'
else:
    log_file = args.logfile

# Load genome digest bed into memory
# N.B. 1 million lines
#bed_col_names = ['chr', 'start', 'stop']
#re_df = pd.read_csv(args.genome_digest_bed, 
#                    header = None, names = bed_col_names, sep = '\t', 
#                    dtype = {'chr':'U','start':'int32','stop':'int32'})
## report chromosomes found
#idx_stats = re_df.groupby('chr').count()

# compare with chr from sam file header

# Define functions
def parse_sam_entry(sam_entry):
    read_name = sam_entry.query_name
    parent_read, pe, slice_number = read_name.split("|")
    ref_name = sam_entry.reference_name
    ref_start = sam_entry.reference_start
    ref_stop = sam_entry.reference_end
    mapped = None
    multimapped = None
    coords = None
    # Check if read mapped
    if sam_entry.is_unmapped:
        mapped = 0
        multimapped = 0
        ref_name = "unmapped"
        ref_start = ""
        ref_stop = ""
        coords = ""
    else:
        mapped = 1
        coords = f'{ref_name}:{ref_start}-{ref_stop}'
        # Check if multimapped
        if sam_entry.is_secondary:
            multimapped = 1
        else:
            multimapped = 0
    return [read_name, parent_read, pe, slice_number, mapped, multimapped, 
            ref_name, ref_start, ref_stop, coords]
    

def process_bam(bam_file):
    # Set up pandas dataframe to store data all sam entries per fragment (read name)
    aln_col_names = ['read_name','parent_read', 'pe', 'slice', 'mapped', 'multimapped',
                         'chr', 'start', 'stop', 'coordinates']
    aln_df = pd.DataFrame(columns = aln_col_names)
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        # loop over the sam file alignment by alignment. Assume sorted by query name
        read_count = 0
        chunk_count = 0
        sam_list = []
        df_list = []
        for aln in bam.fetch(until_eof=True): # until_eof=True required for non co-ordinate sorted and indexed bam
            read_count += 1
            sam_list.append(parse_sam_entry(aln))
            if read_count % chunk_size == 0:
                chunk_count += 1
                print(f'{read_count} sam records processed')
                sam_df = pd.DataFrame(sam_list, columns = aln_col_names)
                df_list.append(sam_df)
                sam_list = []
            # Add current entry to dataframe and update previous read name for next iteration
        sam_df = pd.DataFrame(sam_list, columns = aln_col_names)
        df_list.append(sam_df)
        print(f'Combining {chunk_count} chunks of {chunk_size} sam records')
        aln_df = pd.concat(df_list)
    return aln_df


@numba.jit
def classify_fragments(alignment_df):
    '''Summarise data across all sam entries for a fragment (read_name)'''
    fragment_df = alignment_df.sort_values(['parent_read','coordinates']).groupby(
                'parent_read', as_index=False).agg(
                {'slice':'count', 
                 'pe': 'nunique', 
                 'mapped': 'sum',
                 'multimapped':'sum', 
                 'capture':'nunique', 
                 'capture_counts':'sum',
                 'exclude':'nunique',
                 'exclusion_counts': 'sum',
                 'coordinates':'|'.join})
    fragment_df['capture'] = fragment_df['capture']-1
    fragment_df['exclude'] = fragment_df['exclude']-1
    fragment_df['reporter'] = (fragment_df['mapped'] 
                                - fragment_df['exclusion_counts'] 
                                - fragment_df['capture_counts'])
    return fragment_df


def slice_stats(slice_df):
    stats = slice_df.agg({'read_name': 'nunique',
                         'parent_read': 'nunique',
                         'mapped': 'sum',
                         'multimapped':'sum',
                         'capture': lambda col: col.nunique() - 1,
                         'capture_counts': 'sum',
                         'exclusion_counts': 'sum'})
    return stats

def frag_stats(frag_df):
    stats = frag_df.agg({'parent_read': 'nunique',
                         'mapped': lambda col: (col > 1).sum(),
                         'multimapped': lambda col: (col > 0).sum(),
                         'capture_counts': lambda col: (col > 0).sum(),
                         'exclusion_counts': lambda col: (col > 0).sum(),
                 })
    return stats

def main():
    # Open log file
    logf = open(log_file, 'w')
    # Process bam file
    bam_df = process_bam(args.input_bam)
    bam_time = time.time()
    logf.write(f"Bam file processed: {bam_df.shape[0]} records in {bam_time-start} seconds\n")
    print(f"Bam file processed: {bam_df.shape[0]} records in {bam_time-start} seconds")
    # Intersect with bed files to get capture and exclusion site intersections
    capture_df = pd.read_csv(args.capture, header = None, 
                         names = ["read_name","capture"], sep = '\t', 
                         index_col = 'read_name')
    logf.write(f"Capture file loaded: {capture_df.shape[0]} captured slices\n")
    print(f"Capture file loaded: {capture_df.shape[0]} captured slices")
    cc_df = pd.read_csv(args.capture_count, sep='\t', 
                        header=None, names=['read_name', 'capture_counts'], 
                        index_col = 'read_name')
    logf.write(f"Capture count file loaded: {cc_df.shape[0]} captured slices\n")
    print(f"Capture count file loaded: {cc_df.shape[0]} captured slices")
    exclude_df = pd.read_csv(args.exclude, header = None, 
                         names = ["read_name","exclude"], sep = '\t', 
                         index_col = 'read_name')
    logf.write(f"Exclude file loaded: {exclude_df.shape[0]} excluded slices\n")
    print(f"Exclude file loaded: {exclude_df.shape[0]} excluded slices")
    ex_df = pd.read_csv(args.exclude_count, sep='\t', 
                        header=None, names=['read_name', 'exclusion_counts'], 
                        index_col = 'read_name')
    logf.write(f"Excluded count file loaded: {ex_df.shape[0]} excluded slices\n")
    print(f"Excluded count file loaded: {ex_df.shape[0]} excluded slices")
    bed_time = time.time()
    print(f"All external files loaded in {bed_time-bam_time} seconds")
    # Merge all dataframes
    print(f'Merging dataframes...')
    slice_df = bam_df.merge(capture_df, 
                        on = "read_name", how = "left").merge(exclude_df, 
                        on = "read_name", how = "left").fillna("-")
    slice_df = slice_df.merge(cc_df, 
                        on = "read_name", how = "left").merge(ex_df, 
                        on = "read_name", how = "left").fillna(0)
    # Drop dataframes after merge
    del capture_df, exclude_df, cc_df, ex_df
    merge_time = time.time()
    print(f"Dataframes merged in {merge_time-bed_time} seconds")
    # Write results to file
    write_aln_time = merge_time
    if debug:
        print(f'Writing annotated sam data to file...')
        with open(output_prefix + ".all.slice.tsv", 'w') as aln_out:
            slice_df.to_csv(aln_out, header = True, sep = '\t', index = False)
        write_aln_time = time.time()
        print(f"Dataframe written to file in {write_aln_time-merge_time} seconds")
    # Groupby read name to classify capture-c fragments
    print('Classifying fragments...')
    frag_df = classify_fragments(slice_df)
    classify_time = time.time()
    print(f"Fragments classified in {classify_time-write_aln_time} seconds")
    # Write unfiltered fragment data to file
    write_frag_time = classify_time
    if debug:
        print('Writing fragments to file')
        with open(output_prefix + ".all.frag.tsv", 'w') as frag_out:
            frag_df.to_csv(frag_out, header = True, sep = '\t', index = False)
        write_frag_time = time.time()
        print(f"Fragment dataframe written to file in {write_frag_time-classify_time} seconds")
    # Fragment filtering and stats
    logf.write(f'Unfiltered fragment stats:\n{frag_stats(frag_df)}\n')
    print(f'Unfiltered fragment stats:\n{frag_stats(frag_df)}\n')
    # remove duplicate fragments based on slice coordinates
    print('Removing duplicates')
    frag_df.drop_duplicates(subset="coordinates", 
                                keep = False, inplace = True)
    logf.write(f'After duplicate filtering:\n{frag_stats(frag_df)}\n')
    print(f'After duplicate filtering:\n{frag_stats(frag_df)}\n')
    # Report only fragments with capture site and reporter site
    print('Filtering for useful fragments')
    frag_df.query('capture == 1  and reporter > 0', inplace = True)
    logf.write(f'Useful Fragments:\n{frag_stats(frag_df)}\n')
    print(f'Useful Fragments:\n{frag_stats(frag_df)}\n')
    frag_filter_time = time.time()
    print(f"Fragments filtered in {frag_filter_time-write_frag_time} seconds")
    # Filter slice dataframe to only slices from fragments containing both capture and reporter sites
    logf.write(f'Unfiltered slice stats:\n{slice_stats(slice_df)}\n')
    print(f'Unfiltered slice stats:\n{slice_stats(slice_df)}\n')
    print('Filtering for useful slices')
    filt_slice_df = slice_df[slice_df['parent_read'].isin(frag_df['parent_read'])]
    logf.write(f'After filtering for slices within useful fragments:\n{slice_stats(filt_slice_df)}\n')
    print(f'After filtering for slices within useful fragments:\n{slice_stats(filt_slice_df)}\n')
    # Filter again to remove fragments mapping to excluded regions and unmapped reads
    print('Removing excluded and unmapped fragments')
    filt_slice_df = filt_slice_df.query('exclusion_counts == 0 and mapped == 1')
    logf.write(f'After filtering to remove excluded and unmapped slices\n{slice_stats(filt_slice_df)}\n')
    print(f'After filtering to remove excluded and unmapped slices:\n{slice_stats(filt_slice_df)}\n')
    logf.close()
    slice_filter_time = time.time()
    print(f"Slices filtered in {slice_filter_time-frag_filter_time} seconds")
    # Write useful slices to file
    print(f'Writing useful read slices to file...')
    with open(output_prefix + ".useful.slices.tsv", 'w') as aln_out:
        filt_slice_df.to_csv(aln_out, header = True, sep = '\t', index = False)
    write_aln_time = time.time()
    print(f"Dataframe written to file in {write_aln_time-slice_filter_time} seconds")
    # Filter capture slices and write to bed
    filt_slice_df = filt_slice_df.query('capture_counts == 0')
    slice_bed = filt_slice_df[['chr','start','stop','read_name']]
    print(f'Writing reporter read slices to bed file...')
    with open(output_prefix + ".reporter.bed", 'w') as aln_out:
        slice_bed.to_csv(aln_out, header = False, sep = '\t', index = False)
    write_bed_time = time.time()
    print(f"Dataframe written to bed file in {write_bed_time-write_aln_time} seconds")
    print(f"Total time: {slice_filter_time-start} seconds")

if __name__ == '__main__':
    main()
        