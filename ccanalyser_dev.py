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
p.add_argument('-o', '--output', help='output file output_prefix')
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

def format_dataframe_for_printing(f):
    @wraps(f)
    def wrapped(*args, **kwargs):
        df = f(*args, **kwargs)
        return '\n'.join(str(df).split('\n')[:-1])
    return wrapped

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
                                  })
                             .assign(capture=lambda df: df['capture'] - 1,
                                     exclusion=lambda df: df['exclusion'] -1,
                                     reporter_count=lambda df: df['mapped'] - (df['exclusion_count'] + df['capture_count'] + df['blacklist']))
                             .assign(reporter_count=lambda df: np.where(df['capture_count'] > 0, df['reporter_count'], 0))
                             .rename(columns={'capture': 'unique_capture_sites',
                                              'exclusion': 'unique_exclusion_sites',
                                              'restriction_fragment': 'unique_restriction_fragments',
                                              'slice': 'unique_slices',
                                              'blacklist': 'blacklisted_slices'
                                             }))    
    return df_fragments

@format_dataframe_for_printing
def get_slice_stats(df):
    stats = (df.agg({'read_name': 'nunique',
                     'parent_read': 'nunique',
                     'mapped': 'sum',
                     'multimapped':'sum',
                     'capture': 'nunique',
                     'capture_count': 'sum',
                     'exclusion_count': 'sum',
                     'blacklist': 'sum'})
               .rename({'read_name': 'unique_slices',
                        'parent_read': 'unique_fragments',
                        'multimapped': 'multimapping_slices',
                        'capture': 'unique_capture_sites',
                        'capture_count': 'number_of_capture_slices',
                        'exclusion_count': 'number_of_slices_in_exclusion_region',
                        'blacklist': 'number_of_slices_in_blacklisted_region',
                        })
             )
    return stats

@format_dataframe_for_printing
def get_frag_stats(df):
    stats = (df.agg({'parent_read': 'nunique',
                    'mapped': lambda col: (col > 1).sum(),
                    'multimapped': lambda col: (col > 0).sum(),
                    'capture_count': lambda col: (col > 0).sum(),
                    'exclusion_count': lambda col: (col > 0).sum(),
                    'blacklisted_slices': lambda col: (col > 0).sum(),
                    'reporter_count': lambda col: (col > 0).sum()})
                .rename({'parent_read': 'unique_fragments',
                         'multimapped': 'fragments_with_multimapping_slices',
                         'capture_count': 'fragments_with_capture_sites',
                         'exclusion_count': 'fragments_with_excluded_regions',
                         'blacklisted_slices': 'fragments_with_blacklisted_regions',
                         'reporter_count': 'fragments_with_reporter_slices'})
            )
    return stats

@get_timing(task_name='filtering slices')
def filter_slices(df_fragments, df_slices, output_prefix):

    fragstats_fn = f'{output_prefix}.frag.stats'
    slicestats_fn = f'{output_prefix}.slice.stats'

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
        df_fragments.query('reporter_count > 0', inplace=True)
        fragstats_out.write(f'Valid reporter filtered fragment stats:\n{get_frag_stats(df_fragments)}\n')

        #Useful slice filtering
        print('Filtered for useful slices (removing unmapped, excluded, blacklisted, not paired with a capture and slices mapping to two restriction fragments)')
        df_slices_filt = df_slices[df_slices['parent_read'].isin(df_fragments['parent_read'])]
        df_slices_filt = df_slices_filt.query('mapped == 1 and exclusion_count == 0 and blacklist == 0')
        df_slices_filt = (df_slices_filt.sort_values('capture_count', ascending=False) #Sort by capture sites first
                                        .drop_duplicates(subset=['parent_read', 'restriction_fragment']))  #Stops the same restriction fragment being refered to multiple times
        slicestats_out.write(f'Useful slice stats:\n{get_slice_stats(df_slices_filt.reset_index())}\n')
        
        # Generate .tsv of useful slices
        print(f'Writing useful slices to {output_prefix}.tsv')
        df_slices_filt.to_csv(f'{output_prefix}.useful_slices.tsv', sep='\t')

        # Generate .bed file of capture and reporter slices
        print(f'Writing reporter slices to {output_prefix}.reporter.bed')
        df_slices_capture = (df_slices_filt.reset_index()
                                           .query('capture_count > 0'))
        df_slices_reporter = (df_slices_filt.reset_index()
                                           .query('capture_count == 0'))
        
        df_slices_capture[['chrom', 'start', 'end', 'read_name']].to_csv(f'{output_prefix}.capture.bed', sep='\t', header=None, index=False)
        df_slices_reporter[['chrom', 'start', 'end', 'read_name']].to_csv(f'{output_prefix}.reporter.bed', sep='\t', header=None, index=False)        
        
    return df_fragments, df_slices_capture, df_slices_reporter

def get_reporters_by_capture_site(df_capture, df_reporter):
    

    df_capture_merged = (df_capture[['parent_read', 'capture']]
                                   .merge(df_reporter[['parent_read', 'read_name', 'chrom', 'start', 'end']], on='parent_read'))
                                   
    return {capture_site: reporters for capture_site, reporters in df_capture_merged.groupby('capture')}

@get_timing(task_name='analysis of bam file')
def main():

    output_prefix = get_prefix()

    df_alignment = parse_bam(args.input_bam)
    df_alignment = merge_annotations(df_alignment, args.annotations)
    df_fragments = classify_fragments(df_alignment)

    df_fragments_filt, df_capture, df_reporter = filter_slices(df_fragments, df_alignment, output_prefix)

    # Output all reporter fragments grouped by capture site
    print('Outputing reporter fragments by capture site')
    reporters_by_site = get_reporters_by_capture_site(df_capture, df_reporter)
    for cap_site, reporters in reporters_by_site.items():
        reporters[['chrom', 'start', 'end', 'read_name']].to_csv(f'{output_prefix}_{cap_site}.bed', sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
