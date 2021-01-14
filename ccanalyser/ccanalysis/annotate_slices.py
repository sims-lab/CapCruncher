import os
import sys
import pandas as pd
from pandas import errors
import pybedtools
from pybedtools import BedTool
from joblib import Parallel, delayed
import argparse
import itertools
from multiprocessing import Pool
import numpy as np

from ccanalyser.utils.helpers import is_valid_bed, bed_has_name, bed_has_duplicate_names


def find_intersections(
    a,
    b,
    frac=1e-9,
    column_name="count",
    method="get",
    method_na_values={"get": ".", "count": 0},
):

    if isinstance(a, str):
        a = BedTool(a)
    elif isinstance(a, pd.DataFrame):
        a = BedTool.from_dataframe(a)

    if is_valid_bed(b):

        if method == "get":
            bt_intersections = a.intersect(b, loj=True, f=frac)
        elif method == "count":
            bt_intersections = a.intersect(b, loj=True, c=True, f=frac)
        else:
            raise ValueError("method argument must be in [get|count]")

        return format_intersections(
            intersections=bt_intersections, column_name=column_name
        )

    else:
        return format_intersections(
            intersections=a,
            column_name=column_name,
            failed=True,
            na_value=method_na_values.get(method, np.nan),
        )


def format_intersections(intersections, column_name, failed=False, na_value=0):

    df_intersections = intersections.to_dataframe()

    if not failed:
        # Rename last column to column name
        df_intersections = df_intersections.rename(
            columns={df_intersections.columns[-1]: column_name}
        )

        # Extract only needed columns
        df_intersections = df_intersections[
            ["name", 'chrom', 'start', column_name]
        ]

    else:
        df_intersections[column_name] = na_value
        df_intersections = df_intersections[['name', 'chrom', 'start', column_name]]

    return df_intersections.set_index(["name", 'chrom', 'start'])


def cycle_argument(arg):
    """Allows for the same argument to be stated once but repeated for all files"""
    arg = arg.copy()

    if len(arg) == 1:
        return itertools.cycle(arg)
    else:
        return arg


def main(
    actions=None,
    bed1=None,
    bed2=None,
    colnames=None,
    overlap_fractions=None,
    outfile=None,
    duplicates="remove",
):

    # Verify bed integrity
    assert is_valid_bed(bed1), f"Bed1 - {bed1} is invalid"
    assert bed_has_name(bed1), f"Bed1 - {bed1} does not have a name column"
    # assert bed_has_duplicate_names(bed1), f'Bed1 - {bed1} has duplicates in name column'

    print('Formating  bed file')
    # Make base dataframe
    df_bed = df_bed = BedTool(bed1).to_dataframe()

    print('Dealing with any duplicate names')
    # Deal with duplicates -- currently very harsh
    if duplicates == "remove":
        df_bed = (df_bed.sort_values(["score"], ascending=False)
                        .drop_duplicates(subset="name", keep="first"))
        
        bed1 = BedTool.from_dataframe(
            df_bed[["chrom", "start", "end", "name"]]
        )

    # Run the intersection
    n_actions = len(actions)
    dframes = Parallel(n_jobs=n_actions)(
        delayed(find_intersections)(a=bed1, b=b2, frac=f, column_name=name, method=action)
        for b2, f, name, action in zip(
            bed2,
            cycle_argument(overlap_fractions),
            colnames,
            actions,
        )
    )


    # Merge dataframe with annotations
    df_annotations = (
        df_bed.set_index(['name', 'chrom', 'start'])
        .join(dframes, how="left")
        .reset_index()
        .rename(columns={"name": "slice_name"})
        .drop(columns=['score', 'strand'], errors='ignore')
    )

    # Export to csv
    df_annotations.to_csv(outfile, sep="\t", index=False)


# if __name__ == '__main__':

#     parser = argparse.ArgumentParser()
#     parser.add_argument('-a', '--bed1', nargs='+')
#     parser.add_argument('-b', '--bed2', nargs='+')
#     parser.add_argument('--actions', nargs='+', choices=['get', 'count'])
#     parser.add_argument('-c', '--colnames', nargs='+')
#     parser.add_argument(
#         '-f', '--overlap_fractions', nargs='*', default=1e-9, type=float
#     )
#     parser.add_argument('-o', '--outfile', default='out.tsv.gz')
#     args = parser.parse_args()

#     main(**vars(args))
