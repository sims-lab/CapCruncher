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

def load_tsv(tsv, index_field):

    try:
        df = pd.read_csv(tsv, sep='\t', header=0, index_col=index_field)
        
        if not df.empty:
            return df
    
    except Exception as e:
        print(f'Error with {tsv}: {e}')

def join_tsvs(fnames, index_col, n_processes=8):
    index = format_index_var(index_col)
    dframes = Parallel(n_jobs=n_processes)(delayed(load_tsv)(fn, index) for fn in fnames)
    return dframes[0].join(dframes[1:], how='outer')

def concat_tsvs(fnames, delayed=True):
    
    if is_compressed(fnames):
        df = dd.read_csv(fnames, compression='gzip', blocksize=None, sep='\t')
    else:
        df = dd.read_csv(fnames, sep='\t')
    
    if not delayed:
        return df.compute()
    else:
        return df
  
def main(input_files, output_file, index_field, method='join', n_processes=8):

    if method == 'join':
        df = join_tsvs(input_files, index_col=index_field, n_processes=n_processes)
        df.to_csv(output_file, sep='\t')

    elif method == 'concatenate':
        df = concat_tsvs(input_files, delayed=False)       
        df.to_csv(output_file, sep='\t')
    



    
