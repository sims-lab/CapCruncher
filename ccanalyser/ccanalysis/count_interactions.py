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

def main(slices, outfile=None, only_cis=False, remove_exclusions=False, subsample=False):

    df_slices = pd.read_csv(slices, sep='\t') 

    if only_cis:
        df_slices = df_slices.query('capture_chrom == reporter_chrom')
        print('Removed all non-cis interactions')

    if remove_exclusions:
        df_slices = df_slices.query('(reporter_exclusion == ".") or (capture != reporter_exclusion)')
        print('Removed all excluded regions')
    
    if subsample:
        df_slices = df_slices[df_slices['parent_read'].isin(df_slices['parent_read'].sample(n=subsample))]
        print('Subsampled fragments')
    
    
    print('Grouping at the fragment level')
    fragments = df_slices.groupby('parent_read')


    print('Started counting')
    ligated_rf_counts = count_re_site_combinations(fragments, column='reporter_restriction_fragment')

    with xopen.xopen(outfile, mode='wb', threads=4) as w:
        
        header = '\t'.join(['rf1', 'rf2', 'count']) + '\n'
        w.write(header.encode())

        for (rf1, rf2), count in ligated_rf_counts.items():
            line = '\t'.join([rf1, rf2, str(count)]) + '\n'
            w.write(line.encode())


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-f', '--slices')
#     parser.add_argument('-o', '--outfile', default='out.tsv.gz')
#     args = parser.parse_args()

#     main(**vars(args))   