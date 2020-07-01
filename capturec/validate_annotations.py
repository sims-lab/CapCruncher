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

p = argparse.ArgumentParser()
p.add_argument('-i', '--input_file', help='input .tsv file')
p.add_argument('-o', '--output_file', help='output file name', default='joined.tsv')
p.add_argument('--col_names',
               help='additional col names to expect',
               nargs='+',
               default=[])
p.add_argument('--dtypes',
               help='dtypes of additional columns',
               nargs='+',
               default=[])
args = p.parse_args()

def main():

    expected_columns = ['blacklist', 'capture_count',
                        'capture', 'exclusion_count', 'exclusion',
                        'restriction_fragment', *args.col_names]
    expected_columns_dtypes = ['number', 'number',
                               'object', 'number', 'object', 'object',
                               *args.dtypes]

    df = pd.read_csv(args.input_file,
                     sep='\t',
                     index_col='read_name')

    fill_vals = {'object': '-', 'number': 0}
    for col, dtype in zip(expected_columns, expected_columns_dtypes):
        if not col in df.columns:
            df[col] = fill_vals[dtype]

    df.to_csv(args.output_file, sep='\t')


if __name__ == '__main__':
    main()
