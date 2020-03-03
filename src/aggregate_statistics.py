import glob
import itertools
import os
import sys
import argparse

import numpy as np
import pandas as pd


p = argparse.ArgumentParser()
p.add_argument('-d', '--working_directory')
p.add_argument('--stats')
p.add_argument('--reporters')
args = p.parse_args()


def split_fn(fn_ser):
    return pd.DataFrame({'sample': fn_ser.str.split('.').str[0].str.split('/').str[-1],
                         'read_type': fn_ser.str.split('.').str[1].str.split('_').str[0]
                         })


def aggregate_dedup_stats(fnames):

    dframes = [pd.read_csv(fn,
                           sep='\t',
                           header=None,
                           index_col=0,
                           names=['stat', f'{fn.split("/")[-1].split(".")[0]}']
                           )
               .transpose()
               for fn in fnames]

    df = pd.concat(dframes)

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


def aggregate_digestion_stats(fnames):

    df = pd.concat([pd.read_csv(fn, sep='\t', header=0, index_col=0)
                    .transpose()
                    .assign(sample=fn.split('/')[-1].split('.')[0])
                    for fn in fnames], sort=True)

    df.fillna(0, inplace=True)

    df_summary = (df.reset_index()
                  .rename(columns={'index': 'read_type'})
                  .groupby(['sample', 'read_type'])
                  .sum())

    # Get flashed stats
    flashed = (df_summary['total_read_pairs_processed']
               .drop_duplicates()
               .reset_index()
               .assign(read_type=lambda df: df['read_type'].str.replace('read_1', 'pe'))
               .rename(columns={'total_read_pairs_processed': 'flashed_or_unflashed'}))

    # Calculate reads with restriction sites from histogram
    df_hist = df_summary.iloc[:, :-3]
    df_hist.columns = df_hist.columns.astype(int)

    digested = (df_hist.loc[:, df_hist.columns > 0]
                .sum(axis=1)
                .reset_index()
                .assign(read_type=lambda df: df['read_type'].str.replace('read_1', 'pe')
                        .str.replace('read_2', 'pe'))
                .groupby(['sample', 'read_type'])
                .max()
                .reset_index()
                .rename(columns={0: 'read_pairs_with_restriction_sites'}))

    return (flashed, digested)


def aggregate_ccanalyser_stats(fnames):

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

    df_ccanalyser_stats = pd.concat(dframes)
    ccanalyser_stats = (df_ccanalyser_stats.reset_index()[
                        'unique_fragments'])
    ccanalyser_stats = (df_ccanalyser_stats.reset_index()
                        ['index']
                        .str.split('|', expand=True)
                        .rename(columns={0: 'sample', 1: 'read_type', 2: 'stat_type'}))

    ccanalyser_stats['values'] = ccanalyser_stats
    ccanalyser_stats = (ccanalyser_stats.set_index(['sample', 'read_type'])
                        .pivot(columns='stat_type')
                        ['values']
                        [['mapped', 'contains_single_capture',
                                    'contains_capture_and_reporter', 'duplicate_filtered']]
                        .reset_index())

    return ccanalyser_stats


def aggregate_reporter_stats(fnames):

    df_reporter = pd.concat([pd.read_csv(fn, sep=',', header=0,
                                         names=['capture_probe', 'cis/trans', 'count'])
                             .assign(fn=fn.split('/')[-1])
                             for fn in fnames], sort=True)

    df_reporter = pd.concat([split_fn(df_reporter['fn']),
                             df_reporter], axis=1, ignore_index=True)

    df_reporter_summary = (df_reporter.groupby(['capture_probe', 'cis/trans', 'sample', 'flashed_status'])
                                      .sum()
                                      .reset_index())

    return df_reporter_summary


def combine_stats(stats):

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


def main():

    working_dir = args.working_directory

    dedup_stats_files = glob.glob(f'{working_dir}/deduplicated/*.log')
    digestion_stats_files = glob.glob(f'{working_dir}/digest/*.tsv')
    ccanalyser_stats_files = glob.glob(f'{working_dir}/ccanalyser/stats/*.slice.stats')
    ccanalyser_reporter_stats_files = glob.glob(f'{working_dir}/ccanalyser/stats/*.reporter.stats')

    total_reads, dedup = aggregate_dedup_stats(dedup_stats_files)
    flashed, digested = aggregate_digestion_stats(digestion_stats_files)
    slice_stats = aggregate_ccanalyser_stats(ccanalyser_stats_files)
    reporter_stats = aggregate_reporter_stats(ccanalyser_reporter_stats_files)

    combined_stats = combine_stats([total_reads, dedup, flashed, digested, slice_stats])

    combined_stats.to_csv(args.stats, sep='\t')
    reporter_stats.to_csv(args.reporters, sep='\t')


if __name__ == '__main__':
    main()
