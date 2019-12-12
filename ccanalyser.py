#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 21:11:23 2019

@author: davids

TODO

- don't allow multiple reports for the same RE fragment from the same capture (DS)
    - filter step before output
- output per viewpoint - bedgraph
    - groupby pandas
- stat file output tab delimited (DS)
    - Digestion efficiency - capture + exclusion stats
- report / plots (AS)
    - jupyter notebook
- make trackhub (overlay) (AS)
    - colours in yml

- tri-c support - talk to Marieke
- mnase support
- chunk input fastq files and recombining - big test file
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
p.add_argument('-a', '--annotations', 
               help='Tab-delimited text file containing annotation for each read in bam file')
p.add_argument('-o', '--output', help='output file prefix')
args = p.parse_args()

chunk_size = 100000
debug = False

# commmand line options: -i /ifs/projects/davids/capture-c/pipeline/miseq30/bam/miseq.digest.bam -a /ifs/projects/davids/capture-c/pipeline/miseq30/ccanalyser/annotations.tsv -o test

# assertions - check all input files exist
assert os.path.isfile(args.input_bam), "Input sam file not found"
assert os.path.isfile(args.annotations), "Annotation file  not found"


# Set default output files
output_prefix = None
if args.output is None:
    output_prefix = os.path.basename(args.input_bam)
else:
    output_prefix = args.output

# Set stat file names
slice_stats_file = output_prefix + '.slice.stats'
frag_stats_file = output_prefix + '.frag.stats'


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
                 'capture_count':'sum',
                 'exclusion':'nunique',
                 'exclusion_count': 'sum',
                 'restriction_fragment': 'nunique',
                 'blacklist': 'sum',
                 'coordinates':'|'.join})
    fragment_df['capture'] = fragment_df['capture']-1
    fragment_df['exclusion'] = fragment_df['exclusion']-1
    fragment_df['restriction_fragment'] = fragment_df['restriction_fragment']-1
    fragment_df['reporter'] = (fragment_df['mapped'] 
                                - fragment_df['exclusion_count'] 
                                - fragment_df['capture_count']
                                - fragment_df['blacklist'])
    return fragment_df


def slice_stats(slice_df):
    stats = slice_df.agg({'read_name': 'nunique',
                         'parent_read': 'nunique',
                         'mapped': 'sum',
                         'multimapped':'sum',
                         'capture': 'nunique',
                         'capture_count': 'sum',
                         'exclusion_count': 'sum',
                         'blacklist': 'sum'})
    return stats

def frag_stats(frag_df):
    stats = frag_df.agg({'parent_read': 'nunique',
                         'mapped': lambda col: (col > 1).sum(),
                         'multimapped': lambda col: (col > 0).sum(),
                         'capture_count': lambda col: (col > 0).sum(),
                         'exclusion_count': lambda col: (col > 0).sum(),
                         'blacklist': lambda col: (col > 0).sum(),
                         'reporter': lambda col: (col > 0).sum(),
                 })
    return stats

def main():
    # Process bam file
    bam_df = process_bam(args.input_bam)
    bam_time = time.time()
    print(f"Bam file processed: {bam_df.shape[0]} records in {bam_time-start} seconds")
    # Intersect with bed files to get capture and exclusion site intersections
    annotation_df = pd.read_csv(args.annotations, header = 0, sep = '\t')
    print(f"Annotations file loaded: {annotation_df.shape[0]} annotated slices")
    annotation_time = time.time()
    print(f"Annotations loaded in {annotation_time-bam_time} seconds")
    # Merge all dataframes
    print(f'Merging dataframes...')
    slice_df = bam_df.merge(annotation_df, 
                        on = "read_name", how = "left").fillna("-")
    # Drop dataframe after merge
    del annotation_df
    merge_time = time.time()
    print(f"Dataframes merged in {merge_time-annotation_time} seconds")
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
    # Open frag stats file
    frag_stats = open(frag_stats_file, 'w')
    frag_stats.write(f'Unfiltered fragment stats:\n{frag_stats(frag_df)}\n')
    print(f'Unfiltered fragment stats:\n{frag_stats(frag_df)}\n')
    # remove duplicate fragments based on slice coordinates
    print('Removing duplicates')
    frag_df.drop_duplicates(subset="coordinates", 
                                keep = False, inplace = True)
    frag_stats.write(f'After duplicate filtering:\n{frag_stats(frag_df)}\n')
    print(f'After duplicate filtering:\n{frag_stats(frag_df)}\n')
    # Report only fragments with capture site and reporter site
    print('Filtering for useful fragments')
    frag_df.query('capture == 1  and reporter > 0', inplace = True)
    frag_stats.write(f'Useful Fragments:\n{frag_stats(frag_df)}\n')
    frag_stats.close()
    print(f'Useful Fragments:\n{frag_stats(frag_df)}\n')
    frag_filter_time = time.time()
    print(f"Fragments filtered in {frag_filter_time-write_frag_time} seconds")
    # Filter slice dataframe to only slices from fragments containing both capture and reporter sites
    # Open slice stats file
    slice_stats = open(slice_stats_file, 'w')
    slice_stats.write(f'Unfiltered slice stats:\n{slice_stats(slice_df)}\n')
    print(f'Unfiltered slice stats:\n{slice_stats(slice_df)}\n')
    print('Filtering for useful slices')
    filt_slice_df = slice_df[slice_df['parent_read'].isin(frag_df['parent_read'])]
    slice_stats.write(f'After filtering for slices within useful fragments:\n{slice_stats(filt_slice_df)}\n')
    print(f'After filtering for slices within useful fragments:\n{slice_stats(filt_slice_df)}\n')
    # Filter again to remove fragments mapping to excluded regions and unmapped reads
    print('Removing excluded and unmapped fragments')
    filt_slice_df = filt_slice_df.query('exclusion_count == 0 and mapped == 1')
    slice_stats.write(f'After filtering to remove excluded and unmapped slices\n{slice_stats(filt_slice_df)}\n')
    print(f'After filtering to remove excluded and unmapped slices:\n{slice_stats(filt_slice_df)}\n')
    slice_stats.close()
    slice_filter_time = time.time()
    print(f"Slices filtered in {slice_filter_time-frag_filter_time} seconds")
    # Write useful slices to file
    print(f'Writing useful read slices to file...')
    with open(output_prefix + ".useful.slices.tsv", 'w') as aln_out:
        filt_slice_df.to_csv(aln_out, header = True, sep = '\t', index = False)
    write_aln_time = time.time()
    print(f"Dataframe written to file in {write_aln_time-slice_filter_time} seconds")
    # Filter capture slices and write to bed
    filt_slice_df = filt_slice_df.query('capture_count == 0')
    slice_bed = filt_slice_df[['chr','start','stop','read_name']]
    print(f'Writing reporter read slices to bed file...')
    with open(output_prefix + ".reporter.bed", 'w') as aln_out:
        slice_bed.to_csv(aln_out, header = False, sep = '\t', index = False)
    write_bed_time = time.time()
    print(f"Dataframe written to bed file in {write_bed_time-write_aln_time} seconds")
    print(f"Total time: {slice_filter_time-start} seconds")

if __name__ == '__main__':
    main()
        