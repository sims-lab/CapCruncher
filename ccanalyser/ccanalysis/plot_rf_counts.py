import os
import re
import sys
import argparse

import cooler
import pandas as pd
from pybedtools import BedTool

from ccanalyser.ccanalysis.annotate_slices import bed_has_duplicate_names, bed_has_name



def format_coordinates(coordinates):
    if re.match(r"chr[0-2xXyYmM][0-9]*:\d+-\d+$", coordinates):
        coordinates = f'{coordinates} window_0'
        return BedTool(' '.join(re.split(':|-', coordinates)), from_string=True)

    elif re.match(r'(.*).bed'):
        if bed_has_name(coordinates) and not bed_has_duplicate_names:
            return BedTool(coordinates)
        else:
            return (BedTool(coordinates).to_datamframe()
                                        .assign(name=lambda df: 'window_' + df['index'])
                                        [['chrom', 'start', 'end', 'name']]
                                        .pipe(BedTool.from_dataframe))
    
    else:
        raise ValueError('''Coordinates are not provided in the correct format.
                            Provide coordinates in the form chr[NUMBER]:[START]-END or a .bed file''')
    

def get_bins(bt, bin_size=1000):

    df_bins = (bt.makewindows(w=bin_size, b=bt.fn, i='src')
               .to_dataframe())
    
    df_coords = bt.to_dataframe()

    return df_coords.merge(df_bins, left_on='name', right_on='name', how='inner')

def get_rf_midpoints(cooler_obj, interval):

    coord = f'{interval[0]}:{interval[1]}-{interval[2]}'

    return (cooler_obj.bins()
            .fetch(coord)
            .assign(midpoint_start=lambda df: df['start'] + ((df['end'] - df['start'])  // 2),
                    midpoint_end=lambda df: df['midpoint_start'] + 1)
            .reset_index()
            .rename(columns={'index': 'name'})
            [['chrom', 'midpoint_start', 'midpoint_end', 'name']]
            .pipe(BedTool.from_dataframe)
          )


def main(counts_fname, coordinates, method="tri", bin_size=1000):
    
    counts = cooler.Cooler(counts_fname)
    bt_coordinates = format_coordinates(coordinates=coordinates) 
    df_bins = get_bins(bt=bt_coordinates, bin_size=bin_size)


    print(df_bins)





    # for bin_name, bins in df_bins.groupby('name'):
    #     bins = bins.reset_index(drop=True)





    # for interval in bt_coordinates:
    #     bt_rf_midpoints = get_rf_midpoints(counts, interval)






    # bt_bins = get_bins(bt_coordinates)






    
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--counts_fname', required=True)
    parser.add_argument('-c', '--coordinates', required=True)
    parser.add_argument('--method', default='tri')
    args = parser.parse_args()

    main(**vars(args))
