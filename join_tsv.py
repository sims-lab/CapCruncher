#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 21:49:19 2019

@author: davids

Join any number of tab delimited text files on a single column. 
Index column must have the same header name in each file to be joined.
Performs an outer join on the index column.
Assumes files have headers.
"""

import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('-f', '--index_field', help='shared column name to join on')
p.add_argument('-o', '--output_file', help='output file name', 
               default='joined.tsv')
p.add_argument('-i', '--input_files', nargs = "*", help='at least 2 input files')
args = p.parse_args()

print(args.input_files)
if len(args.input_files) > 1:
    df = pd.read_csv(args.input_files[0], sep='\t', header = 0)
    for list_index in range(1,len(args.input_files)):
        df2join = pd.read_csv(args.input_files[list_index], sep='\t', header = 0)
        df = df.merge(df2join, on = args.index_field,
                        how = "outer")
        df.select_dtypes(include=["int"]).fillna(0, inplace = True)
        df.select_dtypes(include=["object"]).fillna("-", inplace = True)
#write joined df to file
    with open(args.output_file, 'w') as fout:
        df.to_csv(fout, header = True, sep = '\t', index = False)
else:
    print("you must supply at least two input files")