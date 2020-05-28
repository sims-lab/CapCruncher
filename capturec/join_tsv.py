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

import argparse
import pandas as pd
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('-f', '--index_field', help='shared column name to join on')
p.add_argument('-o', '--output_file', help='output file name', 
               default='joined.tsv')
p.add_argument('-i', '--input_files', nargs = "+", help='at least 2 input files')
args = p.parse_args()

def format_index_var(var):
    try:
        return int(var)
    except:
        pass
    return str(var)

        
def replace_na(df):
    numeric_cols = df.select_dtypes(include=["number"]).fillna(0).astype(np.int64)
    df[numeric_cols.columns] = numeric_cols
    string_cols = df.select_dtypes(include=["object"]).fillna("-")
    df[string_cols.columns] = string_cols
    return df


def main():
    
    index_field = format_index_var(args.index_field)
    expected_columns = ['read_name', 'blacklist', 'capture_count', 
                        'capture', 'exclusion_count', 'exclusion', 
                        'restriction_fragment']
    expected_columns_dtypes = ['object', 'number', 'number',
                               'object', 'number', 'object', 'object']
    
    dataframes = []
    for tsv in args.input_files:
        try:
            dataframes.append(pd.read_csv(tsv, 
                                          sep='\t',
                                          header=0,
                                          index_col=index_field))
        except Exception as e:
            print(f'Error with {tsv}: {e}')
            
    
    df_merged = dataframes[0].join(dataframes[1:], how='outer')
    df_merged = replace_na(df_merged)
    
    fill_vals = {'object': '-', 'number': 0}
    for col, dtype in zip(expected_columns, expected_columns_dtypes):
        if not col in df_merged:
            df_merged[col] = fill_vals[dtype]

    df_merged.to_csv(args.output_file, header=True, sep='\t')


if __name__ == '__main__':
    main()