import argparse
import os
import sys
import pandas as pd
from pybedtools import BedTool

p = argparse.ArgumentParser()
p.add_argument('--reporter_tsv', help='File name of concatenated reporter_tsv')
p.add_argument('--re_map', help='Bed file of restriction fragments')
p.add_argument('--output', help='Output file name')
args = p.parse_args()

def main():

    df_reporters = pd.read_csv(args.reporter_tsv, sep='\t')
    df_re_map = pd.read_csv(args.re_map, sep='\t',
                            header=None,
                            names=['re_chrom', 're_start', 're_end', 're_name'])


    df_re_counts = (df_reporters.groupby('restriction_fragment')
                                .size()
                                .reset_index()
                                .rename(columns={0: 'count'})
                    )

    df_re_counts = df_re_counts.merge(df_re_map,
                                      left_on='restriction_fragment',
                                      right_on='re_name')

    bedgraph = (df_re_counts[['re_chrom', 're_start', 're_end', 'count']]
                            .pipe(BedTool.from_dataframe)
                            .sort()
                            .saveas(args.output)
                )


if __name__ == '__main__':
    main()
