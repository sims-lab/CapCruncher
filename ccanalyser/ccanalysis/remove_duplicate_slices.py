import argparse
import os
import sys
import pandas as pd
import dask
from dask.distributed import Client, LocalCluster
import dask.array as da
import dask.dataframe as dd
import random
import mmh3


def mmh3_column(col, hash_type=64):

    if hash_type == 32:
        return [mmh3.hash(x, seed=42) for x in col]
    elif hash_type == 64:
        return [mmh3.hash64(x, seed=42)[0] for x in col]
    elif hash_type == 128:
        return [mmh3.hash128(x, seed=42)[0] for x in col]


def main(
    input_files,
    deduplicated_fragments="deduplicated.hdf5",
    mode="fragments",
    output=None,
    shuffle=False,
    n_cores=8,
    max_memory='64GB',
    sample_name='',
    read_type='',
):

    # client = Client(n_workers=n_cores, 
    #                 threads_per_worker=1,
    #                 memory_limit=max_memory)
    
    

    if mode == "fragments":

        if shuffle:
            random.shuffle(input_files)

        read_csv_options = dict()

        if any(".gz" in fn for fn in input_files):
            read_csv_options["compression"] = "gzip"
            read_csv_options["blocksize"] = None

        deduplicated_ids = (
            dd.read_csv(
                input_files,
                sep="\t",
                usecols=["parent_read", "coordinates"],
                **read_csv_options
            )
            .map_partitions(lambda df: df.apply(mmh3_column)) # Hash to reduce memory
            .drop_duplicates(subset='coordinates')
            .assign(duplicated=0)
            [['parent_read','duplicated']]
            .set_index('parent_read')
            .to_hdf(deduplicated_fragments, key='deduplicated', mode='w'))

    elif mode == "slices":

      

        df_slices = (pd.read_csv(input_files, sep='\t')
                       .assign(parent_read_hashed=lambda df: mmh3_column(df['parent_read']))
                       .set_index('parent_read_hashed'))
        
        n_fragments_total = df_slices['parent_read'].nunique()

        dd_fragments_deduplicated = dd.read_hdf(deduplicated_fragments, key='deduplicated', mode='r')
        df_deduplicated = dd_fragments_deduplicated.join(df_slices, how='inner').compute()
        df_deduplicated.reset_index(drop=True).to_csv(output, sep='\t', index=False)

        n_fragments_unique = df_deduplicated['parent_read'].nunique()

        df_stats = pd.DataFrame()
        df_stats['stat_type'] = ['not-deduplicated', 'deduplicated']
        df_stats['stat'] = [n_fragments_total, n_fragments_unique]
        df_stats['sample'] = sample_name
        df_stats['read_type'] = read_type
        df_stats['stage'] = 'deduplicate_slices'
        

    

    #client.close()


if __name__ == "__main__":

    # cluster = LocalCluster(n_workers=8, threads_per_worker=1)
    # client = Client(cluster)


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
    parser_fragments.add_argument('-p', '--n_cores', default=8, type=int)
    parser_fragments.add_argument('-m', '--max_memory', default='64GB', type=str)

    parser_slices = subparser.add_parser("slices")
    parser_slices.add_argument("-i", "--input_files", required=True)
    parser_slices.add_argument("-f", "--deduplicated_fragments", required=True)
    parser_slices.add_argument("-o", "--output", default="deduplicated.tsv.gz")
    parser_slices.add_argument('-p', '--n_cores', default=8, type=int)
    parser_slices.add_argument('-m', '--max_memory', default='64GB', type=str)

    args = parser.parse_args()

    main(**vars(args))

