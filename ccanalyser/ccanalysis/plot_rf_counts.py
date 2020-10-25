import os
import re
import sys
import argparse
from collections import Counter

import cooler
import pandas as pd
from pybedtools import BedTool
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns

from annotate_slices import bed_has_duplicate_names, bed_has_name


def format_coordinates(coordinates):
    if re.match(r"chr[0-2xXyYmM][0-9]*:\d+-\d+$", coordinates):
        coordinates = f'{coordinates} window_0'
        return BedTool(' '.join(re.split(':|-', coordinates)), from_string=True)

    elif re.match(r'(.*).bed', coordinates):
        if bed_has_name(coordinates) and not bed_has_duplicate_names:
            return BedTool(coordinates)
        else:
            return (BedTool(coordinates).to_dataframe()
                                        .reset_index()
                                        .assign(name=lambda df: 'window_' + df['index'].astype('string'))
                                        [['chrom', 'start', 'end', 'name']]
                                        .pipe(BedTool.from_dataframe))
    
    else:
        raise ValueError('''Coordinates are not provided in the correct format.
                            Provide coordinates in the form chr[NUMBER]:[START]-END or a .bed file''')

def convert_interval_to_coords(interval):
    return f'{interval["chrom"]}:{interval["start"]}-{interval["end"]}'

def get_bins(bt, bin_size=1000):
    return (bt.makewindows(w=bin_size, b=bt.fn, i='src')
               .to_dataframe())

def get_rf_midpoints(cooler_obj, coord):
    return (cooler_obj.bins()
            .fetch(coord)
            .assign(midpoint_start=lambda df: df['start'] + ((df['end'] - df['start'])  // 2),
                    midpoint_end=lambda df: df['midpoint_start'] + 1)
            .reset_index()
            .rename(columns={'index': 'name'})
            [['chrom', 'midpoint_start', 'midpoint_end', 'name']]
            .pipe(BedTool.from_dataframe)
          )

def normalise_counts(df, scale=1000000):
    n_restriction_fragments = df['window_1_n_rf'] * df['window_2_n_rf']
    total = df['total'].values[0]
    counts = df['count'].sum()
    
    return ((counts / n_restriction_fragments) * (scale / total)).values[0]









def main(counts_fname, coordinates, method="tri", bin_size=1000, output_prefix=''):
    
    counts = cooler.Cooler(counts_fname)
    bt_coordinates = format_coordinates(coordinates=coordinates) 
    df_bins = get_bins(bt=bt_coordinates, bin_size=bin_size)
    
    for bin_name, bins in df_bins.groupby('name'):

        # Extract the interval
        interval = bt_coordinates.to_dataframe().loc[lambda df: df['name'] == bin_name].iloc[0]
        coord = convert_interval_to_coords(interval)

        # Format bins dataframe
        bt_bins = (bins.reset_index(drop=True)
                       .reset_index()
                       [['chrom', 'start', 'end', 'index']]
                       .pipe(BedTool.from_dataframe))
        

        # Find rf midpoints in the required region
        bt_rf_midpoints = get_rf_midpoints(counts, coord)

        # Determine which rf (midpoints) overlap the bins
        df_overlaps = (bt_rf_midpoints.intersect(bt_bins, loj=True)
                                      .to_dataframe(disable_auto_names=True, header=None)
                                      .iloc[:, [3,7]]
                                      .rename(columns={3:'rf_id', 7:'bin_id'})
                                      .replace('.', pd.NA)
                                      .dropna())

        # Generate a mapping
        bin_mapping = df_overlaps.set_index('rf_id')['bin_id'].to_dict()

        # Count the number of rf per bin
        n_rf_counts = Counter(bin_mapping.values())

        # Extract counts at the required coordinates
        # Annotate with bin mapping, number of rf, total
        df_rf_binned = (counts.pixels()
                              .fetch(coord)
                              .assign(window_1=lambda df: df['bin1_id'].map(bin_mapping),
                                      window_2=lambda df: df['bin2_id'].map(bin_mapping),
                                      window_1_n_rf=lambda df: df['window_1'].map(n_rf_counts),
                                      window_2_n_rf=lambda df: df['window_2'].map(n_rf_counts),
                                      total=lambda df: df['count'].sum())
            )

        # Normalise count data
        df_normalised = (df_rf_binned.groupby(['window_1', 'window_2'])
                                  .apply(normalise_counts)
                                  .reset_index()
                                  .rename(columns={0:'normalised_count'}))
        df_normalised[['window_1', 'window_2']] = df_normalised[['window_1', 'window_2']].astype(int)

        # Generate a new sparse matrix
        normalised_sparse_matrix =  sparse.coo_matrix((df_normalised['normalised_count'], 
                                                      (df_normalised['window_1'], df_normalised['window_2'])))
        
        # Plot figure
        fig = plt.figure(figsize=(10, 10))
        heatmap = sns.heatmap(data=normalised_sparse_matrix.todense(),
                              vmin=0.01,
                              vmax=20,
                              cmap='viridis',
                              square=True,
                              cbar_kws={"shrink": 0.8},
                              linewidths=0
        )

        # Save figure
        plt.yticks([])
        plt.xticks([])
        #plt.xlabel(coord, fontdict={'size': 12})
        plt.title(f'{interval["name"]} {coord} {bin_size // 1000}kb bins')
        heatmap.figure.savefig(os.path.join(output_prefix, f'{bin_name}.svg'))

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--counts_fname', required=True)
    parser.add_argument('-c', '--coordinates', required=True)
    parser.add_argument('--method', default='tri')
    parser.add_argument('-b', '--bin_size', default=1000, type=int)
    parser.add_argument('-o', '--output_prefix', default='')
    args = parser.parse_args()

    main(**vars(args))
