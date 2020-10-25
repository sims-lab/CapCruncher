import argparse
import os
import sys
import pandas as pd
import dask
import dask.array as da
import dask.dataframe as dd
import random


def save_deduplicate_fragment_ids(
    fragments, output="deduplicated.zarr", key="parent_read"
):

    read_csv_options = dict()

    if any(".gz" in fn for fn in fragments):
        read_csv_options["compression"] = "gzip"
        read_csv_options["blocksize"] = None

    print("Reading fragments and deduplicating")
    deduplicated_ids = (
        dd.read_csv(fragments, sep="\t", **read_csv_options)
        .drop_duplicates(subset="coordinates")
        ["parent_read"]
        .unique()
        .persist()
    )

    print("Calculating maximum id character length")
    max_chars = deduplicated_ids.str.len().max().compute()

    print("Storing ids as zarr array")
    deduplicated_ids = deduplicated_ids.astype(f"S{max_chars}").values
    deduplicated_ids.compute_chunk_sizes()
    deduplicated_ids.to_zarr(output, key, overwrite=True)


def read_deduplicated_fragment_ids(fn, key):
    return da.from_zarr(fn, key).astype(str).compute()


def deduplicate_slices(tsv, dedup_ids):
    df = pd.read_csv(tsv, sep="\t", index_col="parent_read")
    dedup_intersection = df.index.intersection(dedup_ids)
    return df.loc[df.index.isin(dedup_intersection)]


def main(
    input_files,
    deduplicated_fragments="deduplicated.zarr",
    mode="fragments",
    output=None,
    shuffle=False,
):

    if mode == "fragments":

        if shuffle:
            random.shuffle(input_files)

        save_deduplicate_fragment_ids(
            input_files, deduplicated_fragments, key="parent_read"
        )

    elif mode == "slices":
        deduplicated_read_ids = read_deduplicated_fragment_ids(
            deduplicated_fragments, key="parent_read"
        )
        deduplicated_slices = deduplicate_slices(input_files, deduplicated_read_ids)
        deduplicated_slices.to_csv(output, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest="mode")

    parser_fragments = subparser.add_parser("fragments")
    parser_fragments.add_argument("-i", "--input_files", nargs="+", required=True)
    parser_fragments.add_argument("-f", "--deduplicated_fragments", required=True)
    parser_fragments.add_argument(
        "--shuffle",
        help="shuffles the input files to randomise the deduplication",
        action="store_true",
    )

    parser_slices = subparser.add_parser("slices")
    parser_slices.add_argument("-i", "--input_files", required=True)
    parser_slices.add_argument("-f", "--deduplicated_fragments", required=True)
    parser_slices.add_argument("-o", "--output", default="deduplicated.tsv.gz")

    args = parser.parse_args()

    main(**vars(args))
