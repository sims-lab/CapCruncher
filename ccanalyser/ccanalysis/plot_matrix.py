import os
import re
import sys
import argparse
from collections import Counter
from typing import Union

import cooler
from h5py._hl.selections2 import ScalarReadSelection
from iced.filter import filter_high_counts, filter_low_counts
import pandas as pd
from pybedtools import BedTool
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns
import iced
from seaborn.matrix import _matrix_mask

from annotate_slices import bed_has_duplicate_names, bed_has_name
from store_interactions import get_human_readable_number_of_bp, format_restriction_fragment, bin_restriction_fragment

def format_coordinates(coordinates):
    if re.match(r"chr[0-2xXyYmM][0-9]*:\d+-\d+$", coordinates):
        coordinates = f'{coordinates} {coordinates}'
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

def plot_matrix(matrix, figsize=(10, 10), axis_labels=None):
    
    thresh = np.percentile(matrix, 98)
    hm = plt.imshow(matrix, vmax=thresh, vmin=0.0001)
    plt.axis('off')
       
    return hm

def normalise_region_ice(cooler_binned: cooler.Cooler,
                         cooler_restriction_fragments: cooler.Cooler,
                         region: str,
                         filter_low_counts: float = 0,
                         filter_high_counts: float = 0,
                         **normalisation_kwargs):

    matrix = cooler_binned.matrix(balance=False).fetch(region).astype(float)

    if filter_low_counts:
        matrix = iced.filter.filter_low_counts(matrix, percentage=filter_low_counts)
    
    if filter_high_counts:
        matrix = iced.filter.filter_high_counts(matrix)
    


    # Need to remove unwanted keyword args
    del normalisation_kwargs['binning_method']
    del normalisation_kwargs['overlap_fraction']

    return iced.normalization.ICE_normalization(matrix, **normalisation_kwargs)

def extract_bins_from_cooler(cooler_obj, region=None):

    if region:
        return (cooler_obj.bins()
                          .fetch(region)
                          .reset_index()
                          .rename(columns={'index': 'bin'})
                          [['chrom', 'start', 'end', 'bin']])
    else:
        return cooler_obj.bins()[:]

def normalise_region_tric(cooler_binned: cooler.Cooler,
                          cooler_restriction_fragments: cooler.Cooler,
                          region: str,
                          binning_method='overlap',
                          overlap_fraction=0.5,
                          **normalisation_kwargs):

    def normalise_tric_counts(df, scaling_factor=1000000):
        n_restriction_fragments = df['bin1_n_rf'] * df['bin2_n_rf']
        total = df['count'].sum()

        return ((df['count'] / n_restriction_fragments ) * (scaling_factor / total))
    
    bins_binned = extract_bins_from_cooler(cooler_binned, region=region)
    bins_restriction_fragments = extract_bins_from_cooler(cooler_restriction_fragments, region=region)
    bins_restriction_fragments = format_restriction_fragment(bin_restriction_fragment, bin_method=binning_method)

    binned_restriction_fragments = bin_restriction_fragment(bins_binned,
                                                            bins_restriction_fragments,
                                                            overlap_fraction=overlap_fraction)
    
    n_rf_per_bin = (binned_restriction_fragments.groupby('bin')
                                    .size()
                                    .reset_index()
                                    .rename(columns={0: 'rf_per_bin'})
                                    ['rf_per_bin'])
                
    # Determine the number of rf per bin  for the normalisation
    df_sparse_counts = cooler_binned.pixels().fetch(region)
    df_sparse_counts = (df_sparse_counts.merge(n_rf_per_bin, left_on='bin1_id', right_index=True)
                                        .rename(columns={'rf_per_bin': 'bin1_n_rf'})
                                        .merge(n_rf_per_bin, left_on='bin1_id', right_index=True)
                                        .rename(columns={'rf_per_bin': 'bin2_n_rf'})
                                        )
    
    # Normalise count data
    df_sparse_counts['normalised_count'] = normalise_tric_counts(df_sparse_counts)

    # Generate a new sparse matrix
    offset = cooler_binned.offset(region)
    x, y = (df_sparse_counts['bin1_id'] - offset), (df_sparse_counts['bin2_id'] - offset)
    data = df_sparse_counts['normalised_count']
    sm =  sparse.coo_matrix((data, (x, y)))

    return sm.todense()
    

def main(method: str,
         coordinates: str,
         cooler_binned: str,
         cooler_restriction_fragments: str,
         normalisation = None,
         output_prefix: str = None,
         output_format: str = 'png',
         scaling_factor=1e5,
         binning_method='overlap',
         overlap_fraction=0.5,
         filter_low_counts=0.02,
         filter_high_counts=0.02):


    
    cooler_binned = cooler.Cooler(cooler_binned)
    cooler_restriction_fragments = cooler.Cooler(cooler_restriction_fragments)
    resolution =  get_human_readable_number_of_bp(cooler_binned.binsize)

    # Set up possible normalisations
    normalisations = {'ice': normalise_region_ice,
                      'tric': normalise_region_tric}
    
    normalisations_default = {'tri': 'tric',
                              'tiled': 'ice'}
    
    # Select appropriate normalisation
    if normalisation:
        normalisation_func = normalisations[normalisation]
    else:
        normalisation_func = normalisations[normalisations_default[method]]
 
      
    for interval in format_coordinates(coordinates):

        interval_coords = convert_interval_to_coords(interval)
        figure_path = f'{output_prefix}_{resolution}_{interval["name"]}.{output_format}'

        matrix_normalised = normalisation_func(cooler_binned,
                                               cooler_restriction_fragments,
                                               region=interval_coords,
                                               filter_low_counts=filter_low_counts,
                                               filter_high_counts=filter_high_counts,
                                               binning_method=binning_method,
                                               overlap_fraction=overlap_fraction)
        
        fig = plot_matrix(matrix=matrix_normalised)
        fig.figure.savefig(figure_path)
        
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('cooler_restriction_fragments',  help='HDF5 file in Cooler format containing counts of restriction fragment combination counts')
    parser.add_argument('cooler_binned', help='HDF5 file in Cooler format containing binned counts of restriction fragment combination counts')
    parser.add_argument('-c', '--coordinates', required=True, type=str)
    parser.add_argument('--method', choices=['tri', 'tiled'])
    parser.add_argument('--normalisation', default=None)
    parser.add_argument('--scaling_factor', default=1e5, type=int)
    parser.add_argument('--binning_method', default='overlap')
    parser.add_argument('--overlap_fraction', default=0.5, type=float)
    parser.add_argument('--filter_low_counts', default=0.02, type=float)
    parser.add_argument('--filter_high_counts', default=0, type=float)
    parser.add_argument('-o', '--output_prefix', default='')
    parser.add_argument('-f', '--output_format', default='png', choices=['png', 'svg', 'jpeg'])
    args = parser.parse_args()

    main(**vars(args))
