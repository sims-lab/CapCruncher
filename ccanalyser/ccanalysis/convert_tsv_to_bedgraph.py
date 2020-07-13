import argparse
import os
import sys
import pandas as pd
from pybedtools import BedTool

def get_parser(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()
    parser.add_argument('-i','--tsv_input', help='Reporter tsv file')
    parser.add_argument('-b', '--bed', help='''Bed file to intersect with reporters
                                      e.g. RE fragments bed file.''')
    parser.add_argument('--output', help='Output file name', default='out.bedgraph')
    return parser

def main(tsv_input,
         bed,
         output='out.bedgraph'):

    df_reporters = (pd.read_csv(tsv_input, sep='\t') )
    df_reporters[['reporter_start', 'reporter_end']] = df_reporters[['reporter_start', 'reporter_end']].astype(int)
    bt_reporters = BedTool.from_dataframe(df_reporters[['reporter_chrom',
                                                        'reporter_start',
                                                        'reporter_end',
                                                        'reporter_read_name']])

    bt_bed = BedTool(bed)
    bedgraph = (bt_bed.intersect(bt_reporters, c=True)
                      .sort()
                      .cut([0,1,2,4])
                      .saveas(output)
                )

if __name__ == '__main__':
    args = get_parser().parse_args()
    main(**vars(args))
