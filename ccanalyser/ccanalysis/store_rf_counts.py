import os
import argparse
import sys
import pandas as pd
from pybedtools import BedTool
import cooler



def main(rf_counts, rf_map, output='counts.hdf5', only_cis=True):
    # Load counts
    df_rf_counts = pd.read_csv(rf_counts, sep='\t')

    # Load rf map
    df_rf_map = pd.read_csv(rf_map,
                            sep='\t',
                            header=None,
                            names=['chrom', 'start', 'end', 'name'])
    
    # Filter for just rf on the same chromosome as the counts
    if only_cis:
        chromosome = df_rf_counts['rf1'].iloc[0].split('_')[1]
        df_rf_map = df_rf_map.query(f'chrom == "{chromosome}"')
    
    # Generate a mapping of rf name to order in rf dataframe (index)
    # e.g. {chr1_DpnII_0: 0}
    df_rf_map = df_rf_map.reset_index(drop=True)
    rf_mapping = {v: k for k, v in df_rf_map['name'].iteritems()}

    # Generate "pixels" dataframe. Effectively coordinate based sparse matrix
    pixels = (df_rf_counts.assign(bin1_id=df_rf_counts['rf1'].map(rf_mapping).astype(int),
                              bin2_id=df_rf_counts['rf2'].map(rf_mapping).astype(int))
                          .sort_values(['bin1_id', 'bin2_id'])
                          [['bin1_id', 'bin2_id', 'count']])
    
    # Store counts data
    output = output if 'hdf5' in output else f'{output}.hdf5'
    cooler.create_cooler(output,
                         bins=df_rf_map.drop(columns='name'),
                         pixels=pixels[['bin1_id', 'bin2_id', 'count']],
                         ordered=True)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--rf_counts', required=True)
    parser.add_argument('-m', '--rf_map', required=True)
    parser.add_argument('-o', '--output', default='counts.hdf5')
    parser.add_argument('--only_cis', action='store_true', default=False)
    args = parser.parse_args()
    
    main(**vars(args))