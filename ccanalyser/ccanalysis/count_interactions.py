import argparse
import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
from itertools import combinations
import xopen
from natsort import natsort_key

def count_re_site_combinations(fragments, column="restriction_fragment"):

    counts = defaultdict(int)  # Store counts in a default dict

    # For each set of ligated fragments
    for ii, (group_name, frag) in enumerate(fragments):

        if (ii % 10000 == 0) and (ii > 0):
            print(f'Processed {ii} fragments')

        for rf1, rf2 in combinations(frag[column], 2): # Get fragment combinations
            rf1, rf2 = sorted([rf1, rf2], key=natsort_key)     # Sort them to ensure consistency        
            counts[rf1, rf2] += 1

    return counts

def main(slices, outfile=None, remove_exclusions=False, remove_capture=False, subsample=None):

    with xopen.xopen(outfile, mode='wb', threads=4) as writer:
            
        header = '\t'.join(['rf1', 'rf2', 'count']) + '\n'
        writer.write(header.encode())

        df_slices_iterator = pd.read_csv(slices, sep='\t', chunksize=2e6) 

        for ii, df_slices in enumerate(df_slices_iterator):

            if remove_exclusions:
                # Must only remove exclusions if they are relevant for the capture being examined
                print('Removing all excluded regions')
                df_capture = df_slices.query('capture != "."')
                df_reporters_exclusions = df_slices.query('(capture == ".") and (exclusion_count > 0)')

                df_slices_to_remove = (df_capture.drop_duplicates(['parent_read', 'capture'])
                                                    [['parent_read', 'capture']]
                                                    .merge(df_reporters_exclusions[['parent_read', 'exclusion', 'slice_name']], on='parent_read')
                                                    .query('capture == exclusion'))
                
                df_slices = df_slices.loc[~(df_slices['slice_name'].isin(df_slices_to_remove['slice_name']))]
                
        
            if remove_capture:
                df_slices = df_slices.query('capture != "."')

            if subsample:

                if isinstance(subsample, float):
                    subsample_options = {'frac': subsample}
                elif isinstance(subsample, int):
                    subsample_options = {'n': subsample}

                df_slices = df_slices[df_slices['parent_read'].isin(df_slices['parent_read'].sample(**subsample_options))]
                print('Subsampled fragments')
            
            
            print('Grouping at the fragment level')
            fragments = df_slices.groupby('parent_read')


            print('Started counting')
            ligated_rf_counts = count_re_site_combinations(fragments, column='restriction_fragment')


            for (rf1, rf2), count in ligated_rf_counts.items():
                line = '\t'.join([rf1, rf2, str(count)]) + '\n'
                writer.write(line.encode())


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-f', '--slices')
#     parser.add_argument('-o', '--outfile', default='out.tsv.gz')
#     args = parser.parse_args()

#     main(**vars(args))   