import os
import sys
import pandas as pd
import pybedtools
from pybedtools import BedTool
from joblib import Parallel, delayed
import argparse
import itertools
from multiprocessing import Pool


def find_intersections(
    a,
    b,
    frac=1e-9,
    column_name='count',
    method='get',
    method_na_values={'get': '.', 'count': 0},
):

    if is_valid_bed(a):
        
        a = BedTool(a)
        if is_valid_bed(b):
            
            if method == 'get':
                bt_intersections = a.intersect(b, loj=True, f=frac)
            elif method == 'count':
                bt_intersections = a.intersect(b, loj=True, c=True, f=frac)

            return format_intersections(
                bt_intersections, column_name=column_name, bed_named=bed_has_name(a)
            )

        else:
            return format_failed_intersection(
                a, column_name=column_name, na_value=method_na_values[method]
            )

    else:
        raise ValueError(f'{a} is not a valid bed file')


def is_valid_bed(bed):
    try:
        bed = BedTool(bed)
        if bed.field_count(n=1) >= 3:
            return True
    except FileNotFoundError:
        return False


def bed_has_name(bed):
    if isinstance(bed, str):
        bed = BedTool(bed)

    if bed.field_count(n=1) >= 4:
        return True


def format_intersections(intersections, column_name, bed_named=True):

    df_intersections = intersections.to_dataframe()

    if bed_named:
        df_intersections = df_intersections.iloc[:, [3, -1]]
        df_intersections.columns = ['name', column_name]
        return df_intersections.set_index('name')
    else:
        df_intersections = df_intersections.reset_index()
        df_intersections = df_intersections.iloc[:, [0, -1]]
        df_intersections.columns = ['index', column_name]
        return df_intersections.set_index('index')


def format_failed_intersection(bed, column_name, bed_named=True, na_value=0):

    df = bed.to_dataframe()

    if bed_named:
        df = df.iloc[:, 3].to_frame()
        df.columns = ['name']
        df[column_name] = na_value
        return df.set_index('name')
    else:
        df = df.reset_index().iloc[:, 0].to_frame()
        df.columns['index']
        df[column_name] = na_value
        return df.set_index('index')


def format_argument(arg):
    arg = arg.copy()
    
    if len(arg) == 1:
        return itertools.cycle(arg)
    else:
        return arg


def main(
    actions=None,
    bed1=None,
    bed2=None,
    colnames=None,
    overlap_fractions=None,
    outfile=None,
):


    # Uses joblib to perform all operations in parallel
    n_actions = len(actions)
    dframes = Parallel(n_jobs=n_actions)(
        delayed(find_intersections)(a=b1, b=b2, frac=f, column_name=name, method=action)
        for b1, b2, f, name, action in zip(
            format_argument(bed1),
            format_argument(bed2),
            format_argument(overlap_fractions),
            colnames,
            actions,
        )
    )



    if len(bed1) == 1: # Are we just annotating one file
        
        df_annotations = (BedTool(bed1[0])
                          .to_dataframe()
                          ['name']
                         .to_frame()
                         .rename(columns={'name': 'read_name'})
                         .set_index('read_name')
                         )

        df_annotations = df_annotations.join(dframes, how='inner')
        df_annotations = (df_annotations.reset_index()
                                        .rename(columns={'index': 'read_name'})
                        )
        df_annotations.to_csv(outfile, sep='\t', index=False)
    else:
        raise NotImplementedError('Currently just supports one bed1 file')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--bed1', nargs='+')
    parser.add_argument('-b', '--bed2', nargs='+')
    parser.add_argument('--actions', nargs='+', choices=['get', 'count'])
    parser.add_argument('-c', '--colnames', nargs='+')
    parser.add_argument(
        '-f', '--overlap_fractions', nargs='*', default=1e-9, type=float
    )
    parser.add_argument('-o', '--outfile', default='out.tsv.gz')
    args = parser.parse_args()

    main(**vars(args))
