import functools
from loguru import logger
import multiprocessing
import os
import queue
from collections import namedtuple
from multiprocessing import Process
from typing import Iterable, Tuple

import pandas as pd
import ujson
import xxhash
from capcruncher.utils import get_file_type, save_dict


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
                                             Only used if part of a pipeline
         hash_seed (int, optional): Seed to use for hashing. Defaults to 42.
         output_path (os.PathLike, optional): Path to save hashed reads.
        """

        self.inq = inq
        self.hash_seed = hash_seed
        self.output_path = output_path

        super(ReadDeduplicationParserProcess, self).__init__()

    def run(self):
        """Processes fastq reads from multiple files and generates a hashed json dict.

        Dictionary is hashed and in the format {(read  1 name + read 2 name): (s1 + s2)}

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
     hash_seed: Seed for xxhash algorithm. Same as ReadDuplicationParserProcess.
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
         statq (multiprocessing.Queue, optional): Output queue for statistics.
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

        Unique reads are placed on outq and deduplication stats are placed on statq.

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

                        if hash_id not in duplicated_ids:
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


def remove_duplicates_from_parquet(
    slices: Iterable, duplicated_ids: pd.Series, output: os.PathLike
) -> Tuple[int, int]:

    import dask.dataframe as dd
    import pyarrow.dataset as ds

    if not duplicated_ids.empty:
        duplicates = set(duplicated_ids.values)
    else:
        duplicates = set()

    n_reads_total = (
        dd.read_parquet(slices, columns=["parent_id"], engine="pyarrow")["parent_id"]
        .nunique()
        .compute()
    )

    logger.info("Loading and filtering slices")

    # Load and filter data
    slice_dataset = ds.dataset(
        list(slices),
        format="parquet",
    )

    slice_dataset_scanner = slice_dataset.scanner(
        filter=~ds.field("parent_id").isin(duplicates)
    )

    logger.info("Writing unique slices")
    ds.write_dataset(
        slice_dataset_scanner, output, format="parquet", partitioning_flavor="hive"
    )

    n_reads_unique = (
        dd.read_parquet(output, columns=["parent_id"], engine="pyarrow")["parent_id"]
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
