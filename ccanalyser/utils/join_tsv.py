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


def join_tsv(dframes):
    return replace_na(dframes[0].join(dframes[1:], how='outer'))


def concatenate_tsv(dframes):
    return replace_na(pd.concat(dframes))


def main(input_files, output_file, index_field, method='join'):

    index_field = format_index_var(index_field)
    dataframes = []

    for tsv in input_files:
        try:
            df = pd.read_csv(tsv, sep='\t', header=0, index_col=index_field)

            if not df.empty:
                dataframes.append(df)

        except Exception as e:
            print(f'Error with {tsv}: {e}')

    if method == 'join':
        df = join_tsv(dataframes)
        df.to_csv(output_file, header=True, sep='\t')
    elif method == 'concatenate':
        df = concatenate_tsv(dataframes)
        df.to_csv(output_file, header=True, sep='\t')
