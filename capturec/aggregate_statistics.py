import glob
import itertools
import os
import sys
import argparse
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--deduplication_stats', nargs='+',
                   help='Deduplication stats paths' )
    p.add_argument('--digestion_stats', nargs='+',
                   help='Digestion stats paths')
    p.add_argument('--ccanalyser_stats', nargs='+')
    p.add_argument('--reporter_stats', nargs='+')
    p.add_argument('--output_dir')
    return p.parse_args()

def split_fn(fn_ser):
    '''Extracts the sample and read_type attributes from a given file name.

       Args:
        fn_ser: pd.Series containing file names to process

       Returns:
         DataFrame with columns: sample, read_type.
    '''
    return pd.DataFrame({'sample': fn_ser.str.split('.').str[0].str.split('/').str[-1],
                         'read_type': fn_ser.str.split('.').str[1].str.split('_').str[0]
                         })

def combine_dedup_stats(fnames):
    '''Concatenates statistics from initial sequence based read deduplication.

       Args:
        fnames: List/pd.Series containing deduplication stats file names.

      Returns:
        Concatenated dataframe

    '''
    dframes = [pd.read_csv(fn,
                           sep='\t',
                           header=None,
                           index_col=0,
                           names=['stat', f'{fn.split("/")[-1].split(".")[0]}']
                           )
               .transpose()
               for fn in fnames]

    return pd.concat(dframes)


def get_dedup_read_pair_stats(df):
    '''Concatenates statistics from initial sequence based read deduplication.

       Args:
        fnames: List/pd.Series containing deduplication stats file names.

       Returns:
        Concatenated dataframe
    '''

    total_reads = (df['Read_pairs_processed']
                   .reset_index()
                   .assign(read_type='flashed')
                   .rename(columns={'index': 'sample',
                                    'Read_pairs_processed': 'total_reads'})
                   )

    deduplicated = (df['Read_pairs_unique']
                    .reset_index()
                    .assign(read_type='flashed')
                    .rename(columns={'index': 'sample', 'Read_pairs_unique': 'read_pairs_with_unique_sequence'})
                    )

    return (total_reads, deduplicated)


def combine_digestion_stats(fnames):

    dframes = [pd.read_csv(fn).assign(sample=fn.split('/')[-1].split('.')[0])
               for fn in fnames]
    df = pd.concat(dframes)
    return (df.groupby(['bin', 'stat', 'read_type', 'sample'])
              .sum()
              .reset_index())

def get_digestion_read_pair_stats(df):

    # Get flashed stats
    flashed = (df.loc[lambda df: (df['read_type'].isin(['flashed', 'r1'])) &
                                 (df['stat'] == 'total')]
                     .groupby(['sample', 'read_type'])
                     ['frequency']
                     .sum()
                     .reset_index()
                     .assign(read_type=lambda df: df['read_type'].str.replace('r1', 'pe'))
                     .rename(columns={'frequency': 'Flashed or Unflashed'}))


    digested = (df.loc[lambda df: (df['stat'] == 'valid') &
                                  (df['read_type'] != 'r2') &
                                  (df['bin'] != 0)
                      ]
                  .drop(columns='bin')
                  .groupby(['sample', 'read_type'])
                  .sum()
                  .reset_index()
                  [['sample', 'read_type', 'frequency']]
                  .assign(read_type=lambda df: df['read_type'].str.replace('r1', 'pe'))
                  .rename(columns={'frequency': 'read_pairs_with_restriction_site(s)'}))

    return (flashed, digested)


def combine_ccanalyers_stats(fnames):

    fnames_ser = pd.Series(fnames, name='fnames')

    df_fnames = pd.concat([fnames_ser,
                           split_fn(fnames_ser)],
                          axis=1)

    grouped_stats_files = {k: v for k, v in df_fnames.groupby(
        ['sample', 'read_type'])['fnames']}

    dframes = []
    for key, fnames in grouped_stats_files.items():
        df = (pd.concat([pd.read_csv(fn, sep='\t', index_col=0) for fn in fnames])
                .transpose())

        agg_dict = {col: np.sum for col in df.columns.unique()}
        agg_dict.update({'unique_capture_sites': np.max})

        df = pd.concat([pd.Series(agg_dict[k](v, axis=1), name=k)
                        for k, v in df.groupby(df.columns, axis=1)],
                       axis=1)

        df = (df.transpose()
                .add_prefix('|'.join(key) + '|')
                .transpose())

        dframes.append(df)

    return pd.concat(dframes)


def get_ccanalyser_read_pair_stats(df):

    ccanalyser_stats = (df.reset_index()
                        ['index']
                        .str.split('|', expand=True)
                        .rename(columns={0: 'sample', 1: 'read_type', 2: 'stat_type'}))

    ccanalyser_stats['values'] = (df.reset_index()['unique_fragments'])

    ccanalyser_stats = (ccanalyser_stats.set_index(['sample', 'read_type'])
                        .pivot(columns='stat_type')
                        ['values']
                        [['mapped', 'contains_single_capture',
                                    'contains_capture_and_reporter', 'duplicate_filtered']]
                        .reset_index())

    return ccanalyser_stats


def combine_read_pair_stats(stats):

    stats = [df.set_index(['sample', 'read_type']) for df in stats]
    df_stats = stats[0].join(stats[1:], how='outer').fillna(0)

    df_stats_melt = (df_stats.reset_index()
                     .melt(id_vars=['sample', 'read_type'],
                           var_name='stat_type',
                           value_name='read_pairs'))

    df_stats_melt['stat_type'] = (df_stats_melt['stat_type']
                                  .str.capitalize()
                                  .str.replace('_', ' '))

    return df_stats_melt


def combine_reporter_stats(fnames):

    df_reporter = pd.concat([pd.read_csv(fn, sep=',', header=0,
                                         names=['capture_probe', 'cis/trans', 'count'])
                             .assign(fn=fn.split('/')[-1])
                             for fn in fnames], sort=True)



    df_reporter['sample']= df_reporter['fn'].str.split('.').str[0]
    df_reporter['read_type'] = df_reporter['fn'].str.split('.').str[1].str.split('_').str[0]


    return (df_reporter.drop(columns='fn')
                .groupby(['capture_probe', 'cis/trans', 'sample', 'read_type'])
                .sum()
                .reset_index())




def main():


    df_dedup = combine_dedup_stats(args.deduplication_stats)
    df_digestion = combine_digestion_stats(args.digestion_stats)
    df_ccanalyser = combine_ccanalyers_stats(args.ccanalyser_stats)
    df_reporter = combine_reporter_stats(args.reporter_stats)

    stats_dict = dict(zip(['deduplication_stats', 'digestion_stats', 'ccanalyser_stats', 'reporter_stats'],
                          [df_dedup, df_digestion, df_ccanalyser, df_reporter]))


    # Output individual aggregated stats files
    out_dir = args.output_dir.rstrip('/')
    for name, df in stats_dict.items():
        df.to_csv(f'{out_dir}/{name}.tsv', sep='\t')


    # Output combined statst
    total_reads, dedup = get_dedup_read_pair_stats(df_dedup)
    flashed, digested = get_digestion_read_pair_stats(df_digestion)
    slice_stats = get_ccanalyser_read_pair_stats(df_ccanalyser)
    combined_stats = combine_read_pair_stats([total_reads, dedup, flashed, digested, slice_stats])

    combined_stats.to_csv(f'{out_dir}/combined_stats.tsv', sep='\t')


if __name__ == '__main__':
    main(**parse_args())
