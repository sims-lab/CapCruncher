import argparse
import os
import pysam
import pandas as pd
import time
from functools import wraps
import numpy as np
from datetime import timedelta

p = argparse.ArgumentParser()
p.add_argument('-i', '--input_bam', help='BAM file to parse')
p.add_argument('-a', '--annotations', 
               help='Tab-delimited text file containing annotation for each read in bam file')
p.add_argument('-o', '--output', help='output file prefix')
args = p.parse_args()


# assertions - check all input files exist
assert os.path.isfile(args.input_bam), "Input sam file not found"
assert os.path.isfile(args.annotations), "Annotation file  not found"



def get_timing(task_name=None):
    def wrapper(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            time_start = time.perf_counter()
            result = f(*args, **kwargs)
            time_end = time.perf_counter()

            time_taken = timedelta(seconds = (time_end - time_start))
            print(f'Completed {task_name} in {time_taken} (hh:mm:ss.ms)')
            return result
        return wrapped
    return wrapper

def get_prefix():

    if not args.output:
        return os.path.basename(args.input_bam).replace('.bam', '')
    else:
        return args.output

def parse_alignment(aln):
    
    read_name = aln.query_name
    parent_read, pe, slice_number = read_name.split("|")
    ref_name = aln.reference_name
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    # Check if read mapped
    if aln.is_unmapped:
        mapped = 0
        multimapped = 0
        ref_name = "unmapped"
        ref_start = ""
        ref_end = ""
        coords = ""
    else:
        mapped = 1
        coords = f'{ref_name}:{ref_start}-{ref_end}'
        # Check if multimapped
        if aln.is_secondary:
            multimapped = 1
        else:
            multimapped = 0
    return [read_name, parent_read, pe, slice_number, mapped, multimapped, 
            ref_name, ref_start, ref_end, coords]

@get_timing(task_name='processing BAM file')
def parse_bam(bam):
    df_bam = pd.DataFrame([parse_alignment(aln) for aln in pysam.AlignmentFile(bam, 'rb').fetch(until_eof=True)],
                          columns=['read_name','parent_read', 'pe', 'slice', 'mapped', 'multimapped',
                                   'chrom', 'start', 'end', 'coordinates'
                                   ]
                         )
    df_bam.set_index('read_name', inplace=True)
    return df_bam

@get_timing(task_name='merging annotations with BAM input')
def merge_annotations(df, annotations):
    df_ann = pd.read_csv(annotations, sep='\t', header=0, index_col=0)
    return df.join(df_ann, how='inner')

@get_timing(task_name='classifying fragments')
def classify_fragments(df_align):
    
    df_fragments =  (df_align.sort_values(['parent_read','coordinates'])
                             .groupby('parent_read', as_index=False)
                             .agg({'slice':'nunique', 
                                   'pe': 'nunique', 
                                   'mapped': 'sum',
                                   'multimapped':'sum', 
                                   'capture':'nunique', 
                                   'capture_count':'sum',
                                   'exclusion':'nunique',
                                   'exclusion_count': 'sum',
                                   'restriction_fragment': 'nunique',
                                   'blacklist': 'sum',
                                   'coordinates':'|'.join
                                  }
                                 )
                              .assign(capture=lambda df: df['capture'] - 1,
                                     exclusion=lambda df: df['exclusion'] -1,
                                     restriction_fragment=lambda df: df['restriction_fragment'] - 1,
                                     invalid_reporters=lambda df: df['mapped'] - df['exclusion_count'] - df['capture_count'] - df['blacklist']
                                     )
                             .rename(columns={'capture': 'unique_capture_sites',
                                              'exclusion': 'unique_exclusion_sites',
                                              'restriction_fragment': 'unique_restriction_fragments',
                                              'slice': 'unique_slices',
                                              'blacklist': 'blacklisted_slices'
                                             }
                                    )
                  )

    # Valid reporters are defined as mapped reads that are not capture/excusion/blacklist sites
    df_fragments['valid_reporters'] = np.where((df_fragments['capture_count'] != 0) & (df_fragments['mapped'] > 1),
                                              df_fragments['mapped'] - df_fragments['invalid_reporters'],
                                              0)
    
    return df_fragments

def get_slice_stats(df):
    stats =  df.agg({'read_name': 'nunique',
                     'parent_read': 'nunique',
                     'mapped': 'sum',
                     'multimapped':'sum',
                     'capture': 'nunique',
                     'capture_count': 'sum',
                     'exclusion_count': 'sum',
                     'blacklist': 'sum'})
    return stats

def get_frag_stats(df):
    stats = df.agg({'parent_read': 'nunique',
                    'mapped': lambda col: (col > 1).sum(),
                    'multimapped': lambda col: (col > 0).sum(),
                    'capture_count': lambda col: (col > 0).sum(),
                    'exclusion_count': lambda col: (col > 0).sum(),
                    'blacklisted_slices': lambda col: (col > 0).sum(),
                    'valid_reporters': lambda col: (col > 0).sum(),
                 })
    return stats

@get_timing(task_name='filtering slices')
def filter_slices(df_fragments, df_slices):
    
    prefix = get_prefix()
    fragstats_fn = f'{prefix}.frag.stats'
    slicestats_fn = f'{prefix}.slice.stats'

    with open(fragstats_fn, 'w') as fragstats_out,\
         open(slicestats_fn, 'w') as slicestats_out:

        # Unfiltered fragments
        print('Writing unfiltered fragment and slice stats')
        fragstats_out.write(f'Unfiltered fragment stats:\n{get_frag_stats(df_fragments)}\n')
        slicestats_out.write(f'Unfiltered slice stats:\n{get_slice_stats(df_slices.reset_index())}\n')
        

        #Duplicate filtering
        print('Duplicate filtered fragments (looking for exact matches)')
        df_fragments.drop_duplicates(subset="coordinates", keep=False, inplace=True)
        fragstats_out.write(f'Duplicate filtered fragment stats:\n{get_frag_stats(df_fragments)}\n')
        
        #Valid reporter filtering
        print('Filtered for only fragments with valid reporter slices')
        df_fragments.query('valid_reporters > 0', inplace=True)
        fragstats_out.write(f'Valid reporter filtered fragment stats:\n{get_frag_stats(df_fragments)}\n')

        #Useful slice filtering
        print('Filtered for useful slices (removing unmapped, excluded, blacklisted and slices not paired with a capture)')
        df_slices_filt = df_slices[df_slices['parent_read'].isin(df_fragments['parent_read'])]
        df_slices_filt = df_slices_filt.query('mapped == 1 and exclusion_count == 0 and blacklist == 0')
        slicestats_out.write(f'Useful slice stats:\n{get_slice_stats(df_slices.reset_index())}\n')
        
        # Generate .tsv of useful slices
        print(f'Writing useful slices to {prefix}.tsv')
        df_slices_filt.to_csv(f'{prefix}.useful_slices.tsv', sep='\t')

        # Generate .bed file of reporter slices
        print(f'Writing reporter slices to {prefix}.reporter.bed')
        (df_slices_filt.reset_index()
                       .query('capture_count == 0')
                       [['chrom', 'start', 'end', 'read_name']]
                       .to_csv(f'{prefix}.reporter.bed', sep='\t', header=False, index=False)
        )
        
    return df_fragments, df_slices_filt


def main():

    df_alignment = parse_bam(args.input_bam)
    df_alignment = merge_annotations(df_alignment, args.annotations)
    df_fragments = classify_fragments(df_alignment)

    df_fragments_filt, df_alignment_filt = filter_slices(df_fragments, df_alignment)
    





    







if __name__ == '__main__':
    main()
