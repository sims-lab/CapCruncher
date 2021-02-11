import pandas as pd
from collections import defaultdict
from itertools import combinations
import xopen
from ccanalyser.cli import cli
import click

def count_re_site_combinations(fragments, column="restriction_fragment", counts: dict=None):

    if not counts:
        counts = defaultdict(int)  # Store counts in a default dict

    # For each set of ligated fragments
    for ii, (group_name, frag) in enumerate(fragments):

        if (ii % 10000 == 0) and (ii > 0):
            print(f'Processed {ii} fragments')

        for rf1, rf2 in combinations(frag[column], 2): # Get fragment combinations
            rf1, rf2 = sorted([rf1, rf2])     # Sort them to ensure consistency        
            counts[rf1, rf2] += 1

    return counts




@cli.command()
@click.argument('reporters')
@click.option('-o', '--output', help='Name of output file', default='counts.tsv.gz')
@click.option('--remove_exclusions', default=False, help='Prevents analysis of fragments marked as proximity exclusions', is_flag=True)
@click.option('--remove_capture', default=False, help='Prevents analysis of capture fragment interactions', is_flag=True)
@click.option('--subsample', default=0, help='Subsamples reporters before analysis of interactions')
def interactions_count(reporters, output=None, remove_exclusions=False, remove_capture=False, subsample=None):


    '''Counts the number of captured interactions at the restriction fragment level'''

    with xopen.xopen(output, mode='wb', threads=4) as writer:
            
        header = '\t'.join(['bin1_id', 'bin2_id', 'count']) + '\n'
        writer.write(header.encode())

        df_reporters_iterator = pd.read_csv(reporters, sep='\t', chunksize=2e6) 


        ligated_rf_counts = {}
        for ii, df_reporters in enumerate(df_reporters_iterator):

            if remove_exclusions:
                # Must only remove exclusions if they are relevant for the capture being examined
                print('Removing all excluded regions')
                df_capture = df_reporters.query('capture != "."')
                df_reporters_exclusions = df_reporters.query('(capture == ".") and (exclusion_count > 0)')

                df_reporters_to_remove = (df_capture.drop_duplicates(['parent_read', 'capture'])
                                                    [['parent_read', 'capture']]
                                                    .merge(df_reporters_exclusions[['parent_read', 'exclusion', 'slice_name']], on='parent_read')
                                                    .query('capture == exclusion'))
                
                df_reporters = df_reporters.loc[~(df_reporters['slice_name'].isin(df_reporters_to_remove['slice_name']))]
                
        
            if remove_capture:
                df_reporters = df_reporters.query('capture != "."')

            if subsample:

                if isinstance(subsample, float):
                    subsample_options = {'frac': subsample}
                elif isinstance(subsample, int):
                    subsample_options = {'n': subsample}

                df_reporters = df_reporters[df_reporters['parent_read'].isin(df_reporters['parent_read'].sample(**subsample_options))]
                print('Subsampled fragments')
            
            
            print('Grouping at the fragment level')
            fragments = df_reporters.groupby('parent_read')


            print('Started counting')
            ligated_rf_counts = count_re_site_combinations(fragments, column='restriction_fragment', counts=ligated_rf_counts)


            for (rf1, rf2), count in ligated_rf_counts.items():
                line = '\t'.join([str(rf1), str(rf2), str(count)]) + '\n'
                writer.write(line.encode())
