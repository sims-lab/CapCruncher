import argparse
import os
import sys
import pandas as pd
from pybedtools import BedTool

p = argparse.ArgumentParser()
p.add_argument('-i','--input', help='Reporter tsv file')
p.add_argument('-b', '--bed', help='''Bed file to intersect with reporters
                                      e.g. RE fragments bed file.''')
p.add_argument('--output', help='Output file name', default='bedgraph.bedgraph')
args = p.parse_args()

def main():

    df_reporters = (pd.read_csv(args.input, sep='\t') )
    df_reporters[['reporter_start', 'reporter_end']] = df_reporters[['reporter_start', 'reporter_end']].astype(int)
    bt_reporters = BedTool.from_dataframe(df_reporters[['reporter_chrom',
                                                        'reporter_start',
                                                        'reporter_end',
                                                        'reporter_read_name']])

    bt_bed = BedTool(args.bed)


    #Debug
    #bt_reporters.saveas('rep.bed')
    #bt_bed.saveas('bed.bed')

    bedgraph = (bt_bed.intersect(bt_reporters, c=True)
                      .sort()
                      .cut([0,1,2,4])
                      .saveas(args.output)
                )

if __name__ == '__main__':
    main()
