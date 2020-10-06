import argparse
import os
import sys
import pandas as pd
import dask
import dask.dataframe as dd


def get_fragment_coords(df):
    return df.groupby("parent_read").agg(
        coords=pd.NamedAgg(column="reporter_coordinates", aggfunc="|".join)
    )


def main(input_files, output="test.tsv", low_memory=False):

    df = dd.read_csv(
        input_files, compression="gzip", blocksize=None, sep="\t"
    ).persist()

    fragment_coords = (
        df.map_partitions(get_fragment_coords)
        .drop_duplicates(subset=["coords"])
        .compute()
    )

    df.loc[df["parent_read"].isin(fragment_coords.index)].to_csv(
        output, single_file=True, sep="\t"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_files", nargs='+')
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    main(**vars(args))
