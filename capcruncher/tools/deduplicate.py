import functools
import logging
import multiprocessing
import os
import pickle
import queue
from collections import namedtuple
from multiprocessing import Process, Queue, SimpleQueue
from typing import Iterable, List, Literal, NamedTuple, Tuple, Union

import pandas as pd
import ujson
import xxhash
from capcruncher.utils import get_file_type, hash_column, save_dict
from xopen import xopen


class ReadDeduplicationParserProcess(Process):
    """
    Process subclass for parsing fastq file(s) into a hashed {id:sequence} json format.

    Attributes:
     inq: Input read queue
     outq: Output read queue (Not currently used)
     hash_seed: Seed for xxhash64 algorithm to ensure consistency
     save_hash_dict_path: Path to save hashed dictionary
    """

    def __init__(
        self,
        inq: multiprocessing.Queue,
        hash_seed: int = 42,
        output_path: os.PathLike = "parsed.json",
    ):
        """
        Args:
         inq (multiprocessing.SimpleQueue): Input queue for fastq reads.
         outq (multiprocessing.SimpleQueue): Output queue for processed reads.
                                             Only used if part of a pipeline of reader -> parser -> writer
         hash_seed (int, optional): Seed to use for hashing. Defaults to 42.
         save_hashed_dict_path (os.PathLike, optional): Path to save hashed reads in json format. Defaults to "parsed.json".
        """

        self.inq = inq
        self.hash_seed = hash_seed
        self.output_path = output_path

        super(ReadDeduplicationParserProcess, self).__init__()

    def run(self):
        """Processes fastq reads from multiple files and generates a hashed json dictionary.

        Dictionary is hashed and in the format {(read  1 name + read 2 name): (sequence 1 + sequence 2)}.

        Output path is specified by save_hashed_dict_path.

        """

        hash_seed = self.hash_seed
        hash_function = functools.partial(xxhash.xxh64_intdigest, seed=hash_seed)
        records = dict()

        while True:

            try:
                reads = self.inq.get(block=True, timeout=0.01)

                if reads:

                    for read_set in reads:
                        hash_sequence = hash_function(
                            "".join([r.sequence for r in read_set])
                        )
                        hash_id = hash_function("".join([r.name for r in read_set]))
                        records[hash_id] = hash_sequence

                else:
                    break

            except queue.Empty:
                continue
        
        output_format = get_file_type(self.output_path)
        save_dict(records, self.output_path, output_format)


RemovalStatistics = namedtuple(
    "RemovalStatistics", ["reads_total", "reads_unique", "reads_removed"]
)


class ReadDuplicateRemovalProcess(Process):
    """
    Process subclass for parsing fastq file(s) and removing identified duplicates.

    Attributes:
     inq: Input read queue
     outq: Output queue for deduplicated reads.
     duplicated_ids: Concatenated read ids to remove from input fastq files.
     statq: Output queue for statistics.
     reads_total: Number of fastq reads processed.
     reads_unique: Number of non-duplicated reads output.
     hash_seed: Seed for xxhash algorithm. MUST be the same as used by ReadDuplicationParserProcess.
    """

    def __init__(
        self,
        inq: multiprocessing.Queue,
        outq: multiprocessing.Queue,
        stats_tx: multiprocessing.Pipe,
        duplicated_ids: set,
        hash_seed: int = 42,
        hash_read_name: bool = True,
    ):
        """
        Args:
         inq (multiprocessing.SimpleQueue): Input queue for reads to be deduplicated.
         outq (multiprocessing.SimpleQueue): Output queue for deduplicated reads.
         duplicated_ids (set): Hashed read ids to be removed if encountered.
         statq (multiprocessing.Queue, optional): Output queue for statistics. Defaults to None.
         hash_seed (int, optional): Seed for xxhash algorithm. Defaults to 42.
        """

        self.inq = inq
        self.outq = outq
        self.hash_seed = hash_seed
        self.duplicated_ids = duplicated_ids

        # Misc
        self.hash_read_name = hash_read_name

        # Stats
        self.stats_tx = stats_tx
        self.reads_total = 0
        self.reads_unique = 0

        super(ReadDuplicateRemovalProcess, self).__init__()

    def run(self):

        """Performs read deduplication based on sequence.

        Unique reads are placed on outq and deduplication statistics are placed on statq.

        """

        hash_seed = self.hash_seed
        hash_read_name = self.hash_read_name
        hash_function = functools.partial(xxhash.xxh64_intdigest, seed=hash_seed)
        duplicated_ids = self.duplicated_ids
        reads_unique = list()

        while True:

            try:
                reads = self.inq.get(block=True, timeout=0.01)

                if reads:
                    for read_glob in reads:

                        hash_id = hash_function("".join([r.name for r in read_glob]))

                        if not hash_id in duplicated_ids:
                            if hash_read_name:
                                for r in read_glob:
                                    r.name = str(hash_function(r.name))

                            reads_unique.append(read_glob)

                    self.reads_total += len(reads)
                    self.reads_unique += len(reads_unique)
                    self.outq.put(reads_unique.copy())
                    reads_unique.clear()

                else:
                    break

            except queue.Empty:
                continue

        stats = RemovalStatistics(
            self.reads_total, self.reads_unique, self.reads_total - self.reads_unique
        )
        self.stats_tx.send(stats)


def extract_fragment_coordinates(df):

    try:
        coords = df["coordinates"].str.extract(
            r"^chr(?P<chrom1>.*?):(?P<start>\d+).*\|chr(?P<chrom2>.*?):\d+-(?P<end>\d+)"
        )
        coords = coords["chrom1"].str.cat(coords.iloc[:, 1:])
        coords.name = "coordinates_pe"
    except AttributeError as e:  # Might not be any data
        coords = pd.Series(data=[], name="coordinates_pe", dtype="object")

    return coords


def identify_coordinate_duplicates_from_tsv(
    fragments: pd.DataFrame, read_type: Literal["flashed", "pe"], buffer: float = 1e6
):

    df_fragments = pd.read_csv(
        fragments,
        sep="\t",
        chunksize=buffer,
        usecols=["parent_read", "coordinates"],
    )

    unique_coordinates = set()
    duplicated_fragments = set()

    for ii, df in enumerate(fragments):

        # Shuffle to stop fragments at the end of the sample always being removed
        df = df.sample(frac=1)

        if read_type == "flashed":

            lst_id = hash_column(df["parent_read"])
            lst_coords = hash_column(df["coordinates"])

        else:
            # Extract chrom1 + start1 + chrom(last entry) + end(last entry)
            coords_df = df["coordinates"].str.extract(
                r"^chr(?P<chrom1>.*?):(?P<start>\d+).*\|chr(?P<chrom2>.*?):\d+-(?P<end>\d+)"
            )

            lst_id = hash_column(df["parent_read"])

            # chrom1+start1+chrom-1+end-1(hashed)
            lst_coords = hash_column(coords_df["chrom1"].str.cat(coords_df.iloc[:, 1:]))

        for id, coord in zip(lst_id, lst_coords):
            if not coord in unique_coordinates:
                unique_coordinates.add(coord)
            else:
                duplicated_fragments.add(id)

    return duplicated_fragments


def identify_coordinate_duplicates_from_hdf5(
    fragments: list, read_type: Literal["flashed", "pe"]
):

    import dask.dataframe as dd

    try:
        df_fragments_coords = dd.read_hdf(
            fragments, key=f"/fragments", columns=["id", "coordinates"]
        )

        if read_type == "flashed":

            duplicated_ids = (
                df_fragments_coords.map_partitions(
                    lambda df: df.assign(coordinates=hash_column(df["coordinates"]))
                )
                .shuffle(on="coordinates")
                .map_partitions(lambda df: df[df.duplicated(subset="coordinates")])[
                    "id"
                ]
            )
        elif read_type == "pe":

            duplicated_ids = (
                df_fragments_coords.map_partitions(
                    lambda df: df.join(extract_fragment_coordinates(df))[
                        ["id", "coordinates_pe"]
                    ]
                    # .assign(coordinates_pe=lambda df: hash_column(df["coordinates_pe"]))
                )
                .shuffle(on="coordinates_pe")
                .map_partitions(lambda df: df[df.duplicated(subset="coordinates_pe")])[
                    "id"
                ]
            )

    except KeyError as e:
        logging.warn("{e}")
        duplicated_ids = dd.from_pandas(pd.Series(data=[], name="id"), npartitions=1)

    return duplicated_ids


def identify_coordinate_duplicates_from_parquet(
    fragments: list, read_type: Literal["flashed", "pe"]
):

    import dask.dataframe as dd

    df_fragments_coords = dd.read_parquet(
        fragments,
        columns=["id", "coordinates"],
        chunksize="100MB",
        aggregate_files=True,
        engine="pyarrow",
    )

    if read_type == "flashed":

        ids_duplicated = df_fragments_coords.shuffle(on="coordinates").map_partitions(
            lambda df: df[df.duplicated(subset="coordinates")]
        )["id"]

    elif read_type == "pe":

        ids_duplicated = (
            df_fragments_coords.map_partitions(
                lambda df: df.join(extract_fragment_coordinates(df))[
                    ["id", "coordinates_pe"]
                ]
            )
            .shuffle(on="coordinates_pe")
            .map_partitions(lambda df: df[df.duplicated(subset="coordinates_pe")])["id"]
        )

    return ids_duplicated


def remove_duplicates_from_tsv(
    slices: os.PathLike, output: os.PathLike, duplicated_ids: os.PathLike, buffer=1e6
):

    # Remove output if it exist as will need to append to file.
    if os.path.exists(output):
        os.unlink(output)

    df_slices = pd.read_csv(slices, sep="\t", chunksize=buffer)
    ids_duplicated = read_duplicated_ids(duplicated_ids)

    n_reads_total = 0
    n_reads_unique = 0

    # Iterate slices in chunks
    for ii, df in enumerate(df_slices):

        logging.info(f"Processed {(ii + 1) * buffer} slices")

        n_reads_total += df["parent_read"].nunique()

        # Hash the parent_read column and remove any duplicated ids.
        df = (
            df.assign(parent_read_hashed=lambda df: hash_column(df["parent_read"]))
            .set_index("parent_read_hashed")
            .loc[lambda df: ~df.index.isin(ids_duplicated)]
        )

        n_reads_unique += df["parent_read"].nunique()

        # Append to file.
        df.reset_index(drop=True).to_csv(
            output, sep="\t", index=None, header=True if ii < 1 else False, mode="a"
        )

    return (n_reads_total, n_reads_unique)


def remove_duplicates_from_hdf5(
    slices: Iterable, duplicated_ids: pd.Series, output: os.PathLike
) -> Tuple[int, int]:

    import dask.dataframe as dd

    n_slices_total = 0

    try:
        # Need to get total number of slices
        for slice_file in slices:
            with pd.HDFStore(slice_file, "r") as store:
                n_slices_total += store.get_storer("slices").nrows

        ddf = dd.read_hdf(slices, "slices").map_partitions(
            lambda df: df.loc[~(df["parent_id"].isin(duplicated_ids))]
        )

        ddf.to_hdf(
            output,
            key="slices",
            format="table",
            data_columns=["viewpoint"],
            mode="w",
            min_itemsize={
                "slice_name": 75,
                "parent_read": 75,
                "coordinates": 75,
                "chrom": 25,
            },
            complib="blosc",
            complevel=2,
        )

        # Need to get final number of slices
        with pd.HDFStore(output, "r") as store:
            n_slices_unique = store.get_storer(f"slices").nrows

    except KeyError as e:
        # Obviously missing any data (due to filtering)
        n_slices_total = 0
        n_slices_unique = 0
    return (n_slices_total, n_slices_unique)


def remove_duplicates_from_parquet(
    slices: Iterable, duplicated_ids: pd.Series, output: os.PathLike
) -> Tuple[int, int]:

    import dask.dataframe as dd

    duplicates = tuple(duplicated_ids.values)

    n_reads_total = (
        dd.read_parquet(slices, columns=["parent_id"], engine="pyarrow")
        ["parent_id"]
        .nunique()
        .compute()
    )

    logging.info("Loading and filtering slices")
    # Load and filter data
    ddf = dd.read_parquet(
        slices,
        filters=[("parent_id", "not in", duplicates)],
        engine="pyarrow-dataset",
        
    )

    logging.info("Writing unique slices")
    ddf.to_parquet(output, compression="snappy", engine="pyarrow")

    n_reads_unique = (
        dd.read_parquet(output, columns=["parent_id"], engine="pyarrow")
        ["parent_id"]
        .nunique()
        .compute()
    )
    return (n_reads_total, n_reads_unique)


def read_duplicated_ids(path: os.PathLike):

    from xopen import xopen

    file_type = get_file_type(path)

    if file_type == "json":
        with xopen.xopen(path, "r") as r:
            ids_duplicated = {int(x) for x in ujson.load(r)}

    elif file_type == "hdf5":

        try:
            ids_duplicated = pd.read_hdf(path, key="/duplicated_ids")
        except KeyError:
            ids_duplicated = pd.Series(data=["NO_DATA"], name="/duplicated_ids")

    elif file_type == "pickle":
        ids_duplicated = pd.read_pickle(path)

    return ids_duplicated
