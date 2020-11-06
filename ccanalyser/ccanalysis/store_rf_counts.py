import os
import argparse
import sys
import pandas as pd
from pybedtools import BedTool
import cooler


def get_bins(genome='mm9', bin_size=1000) -> pd.DataFrame:

    #TODO: Allow custom genome to be used
    return (BedTool().makewindows(w=bin_size, genome=genome, i='srcwinnum')
               .to_dataframe()
               .reset_index()
               .rename(columns={'index':'bin'}))


def store_counts_fragments(df_counts, df_map, outfile='rf_counts.hdf5'):
    
    # Generate a mapping of rf name to order in rf dataframe (index)
    # e.g. {chr1_DpnII_0: 0}
    df_rf_map = df_map.reset_index(drop=True)
    rf_mapping = {v: k for k, v in df_map['name'].iteritems()}

    # Generate "pixels" dataframe. Effectively coordinate based sparse matrix
    pixels = (df_counts.assign(bin1_id=df_counts['rf1'].map(rf_mapping).astype(int),
                               bin2_id=df_counts['rf2'].map(rf_mapping).astype(int))
                          .sort_values(['bin1_id', 'bin2_id'])
                          [['bin1_id', 'bin2_id', 'count']])
    
    cooler.create_cooler(outfile,
                         bins=df_rf_map.drop(columns='name'),
                         pixels=pixels[['bin1_id', 'bin2_id', 'count']],
                         ordered=True)


def store_counts_binned(df_counts, df_map, genome, bin_size=1000, outfile='binned_counts.hdf5'):
    
    # Get bins
    df_bins = get_bins(genome=genome, bin_size=bin_size)

    # Generate BedTool objects for the bins and rf map
    bt_bins = BedTool.from_dataframe(df_bins[['chrom', 'start', 'end', 'bin']])
    bt_map = BedTool.from_dataframe(df_map[['chrom', 'start', 'end', 'name']])

    # Intersect bins with restriction fragments
    df_bins_rf = (bt_bins.intersect(bt_map, loj=True)
                         .to_dataframe()
                         [['name', 'thickEnd']]
                         .rename(columns={'name': 'bin', 'thickEnd': 'rf'}))
    
    # Merge the new window id with the counts dataframe
    df_counts = (df_counts.merge(df_bins_rf, left_on='rf1', right_on='rf')
                          .merge(df_bins_rf, left_on='rf2', right_on='rf')
                          .rename(columns={'bin_x': 'bin1_id', 'bin_y': 'bin2_id'})
                          .sort_values(['bin1_id', 'bin2_id'])
                          .groupby(['bin1_id', 'bin2_id'])
                          ['count']
                          .sum()
                          .reset_index())


    cooler.create_cooler(outfile,
                         bins=df_bins[['chrom', 'start', 'end']],
                         pixels=df_counts[['bin1_id', 'bin2_id', 'count']],
                         ordered=True)
    

def main(rf_counts, rf_map, genome, binsize=1000, output='counts.hdf5', only_cis=True):
    
    # Load counts
    df_rf_counts = pd.read_csv(rf_counts, sep='\t')

    # Load rf map
    df_rf_map = pd.read_csv(rf_map,
                            sep='\t',
                            header=None,
                            names=['chrom', 'start', 'end', 'name'])
    
    # Filter for just rf on the same chromosome as the counts
    # if only_cis:
    #     chromosome = df_rf_counts['rf1'].iloc[0].split('_')[1]
    #     df_rf_map = df_rf_map.query(f'chrom == "{chromosome}"')
    

    # Store counts data
    output = output if 'hdf5' in output else f'{output}.hdf5'
    store_counts_fragments(df_counts=df_rf_counts, df_map=df_rf_map, outfile=output.replace('.hdf5', '_rf.hdf5'))
    store_counts_binned(df_counts=df_rf_counts,
                       df_map=df_rf_map,
                       genome=genome,
                       bin_size=binsize, 
                       outfile=output.replace('.hdf5', f'_{binsize // 1000}kb_binned.hdf5'))




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--rf_counts', required=True)
    parser.add_argument('-m', '--rf_map', required=True)
    parser.add_argument('-g', '--genome', default='mm9', required=True)
    parser.add_argument('-b', '--binsize', default=1000, type=int)
    parser.add_argument('-o', '--output', default='counts.hdf5')
    parser.add_argument('--only_cis', action='store_true', default=False)
    args = parser.parse_args()
    
    main(**vars(args))