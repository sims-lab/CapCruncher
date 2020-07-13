import os
import sys
import pandas as pd
import argparse

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 21:49:19 2019

@author: asmith

Script validates that all columns that should be present in annotations for
ccanalyser.py exist in the dataframe. This prevents dataframe missing key errors.
Additional columns can be added.

"""
def get_parser(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='input .tsv file')
    parser.add_argument('-o', '--output_file', help='output file name', default='joined.tsv')
    parser.add_argument('--col_names',
                   help='additional col names to expect',
                   nargs='+',
                   default=[])
    parser.add_argument('--dtypes',
                   help='dtypes of additional columns',
                   nargs='+',
                   default=[])

    return parser

def main(input_file,
         output_file='out.tsv.gz',
         col_names=[],
         dtypes=[]):

    expected_columns = ['blacklist', 'capture_count',
                        'capture', 'exclusion_count', 'exclusion',
                        'restriction_fragment', *col_names]
    expected_columns_dtypes = ['number', 'number',
                               'object', 'number', 'object', 'object',
                               *dtypes]

    df = pd.read_csv(input_file,
                     sep='\t',
                     index_col='read_name')

    fill_vals = {'object': '-', 'number': 0}
    for col, dtype in zip(expected_columns, expected_columns_dtypes):
        if not col in df.columns:
            df[col] = fill_vals[dtype]

    df.to_csv(output_file, sep='\t')


if __name__ == '__main__':

    args = get_parser().parse_args()
    main(**vars(args))
