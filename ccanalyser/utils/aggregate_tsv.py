#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 21:49:19 2019

@author: davids, asmith

Join any number of tab delimited text files on a single column.
Index column must have the same header name in each file to be joined.
Performs an outer join on the index column.
Assumes files have headers.
"""
import dask
import argparse
from dask.dataframe.io.csv import to_csv
import numpy as np
import dask.array as da
import dask.dataframe as dd
import pandas as pd
from joblib import Parallel, delayed

def format_index_var(var):

    if var:
        try:
            return int(var)
        except:
            return str(var)
    else:
        return None

def replace_na(df):
    numeric_cols = df.select_dtypes(include=["number"]).fillna(0).astype(np.int64)
    df[numeric_cols.columns] = numeric_cols
    string_cols = df.select_dtypes(include=["object"]).fillna(".")
    df[string_cols.columns] = string_cols
    return df

def is_compressed(files):
    if any('.gz' in fn for fn in files):
        return True

def load_tsv(tsv, index=None, header=None):

    try:
        df = pd.read_csv(tsv, sep='\t', header=header, index_col=index)
        
        if not df.empty:
            return df
    
    except Exception as e:
        print(f'Error with {tsv}: {e}')

def join_tsvs(fnames, index_col, n_processes=8, header=True):
    index = format_index_var(index_col)

    dframes = Parallel(n_jobs=n_processes)(delayed(load_tsv)(fn, index=index, header=header) for fn in fnames)
    return dframes[0].join(dframes[1:], how='outer')

def concat_tsvs(fnames, delayed=False, header=None):

    if is_compressed(fnames):
        df = dd.read_csv(fnames, compression='gzip', blocksize=None, sep='\t', header=header)
    else:
        df = dd.read_csv(fnames, sep='\t', header=header)
    
    return df
    

  
def main(input_files,
         output,
         index=None,
         header=None,
         method=None,
         n_processes=8,
         groupby_columns=None,
         aggregate_method=None,
         aggregate_columns=None
         ):

    header = 0 if header else None

    if method == 'join':
        df = join_tsvs(input_files, index_col=index, n_processes=n_processes)
        df.to_csv(output, sep='\t')

    elif method == 'concatenate':
        df = concat_tsvs(input_files, header=header, delayed=False)
        df.to_csv(output,
                  sep='\t',
                  index=False,
                  single_file=True,
                  compression='gzip' if output.endswith('.gz') else None)
    
    elif method == 'aggregate':
        
        df = load_tsv(input_files, index=index, header=header)

        if groupby_columns:

            agg_dict = dict(zip(aggregate_columns, aggregate_method))

            (df.groupby(groupby_columns)
               .agg(agg_dict)
               .to_csv(output, sep='\t')
            )
        
        else:
            raise ValueError('Groupby columns not provided!')



#if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # subparser = parser.add_subparsers(dest="method")

    # parser_join = subparser.add_parser("join")
    # parser_join.add_argument("-i", "--input_files", nargs="+", required=True)
    # parser_join.add_argument("--index", default='parent_read', required=True)
    # parser_join.add_argument("--header", default=False, action='store_true')
    # parser_join.add_argument('-o',"--output", default='joined.tsv.gz')
    # parser_join.add_argument('-p',"--n_processes", default=8, type=int)

    # parser_concatenate = subparser.add_parser("concatenate")
    # parser_concatenate.add_argument("-i", "--input_files", nargs="+", required=True)
    # #parser_concatenate.add_argument("--index", default=None)
    # parser_concatenate.add_argument("--header", default=False, action='store_true')
    # parser_concatenate.add_argument('-o',"--output", default='concatenated.tsv.gz')

    # parser_aggregate = subparser.add_parser("aggregate")
    # parser_aggregate.add_argument("-i", "--input_files", required=True)
    # parser_aggregate.add_argument("--index", default=None)
    # parser_aggregate.add_argument("--header", default=False, action='store_true')
    # parser_aggregate.add_argument('-o',"--output", default='aggregated.tsv.gz')
    # parser_aggregate.add_argument('-g',"--groupby_columns", nargs='+')
    # parser_aggregate.add_argument('--aggregate_method', nargs='+')
    # parser_aggregate.add_argument('--aggregate_columns', nargs='+')

    # args = parser.parse_args()

    # main(**vars(args))



    
