#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""
import pickle
import sys
from multiprocessing import Process, Queue, SimpleQueue
from typing import Union
import pandas as pd
from collections import Counter
import mmh3
from ccanalyser.utils.io import FastqReaderProcess, FastqWriterProcess

def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return open(fn, "w")
    else:
        return fn


def merge_dictionaries(dicts: list):
    dict_merged = dict()
    for d in dicts:
        dict_merged.update(d)
    return dict_merged


def invert_dictionary(d):
    return {v: k for k, v in d.items()}


class DeduplicationStatistics():
    def __init__(self,
                 sample,
                 read_type: str = 'pe',
                 reads_total=0,
                 reads_unique=0):

        self.sample = sample
        self.read_type = read_type
        self.reads_total = reads_total
        self.reads_unique = reads_unique
        self.reads_removed = reads_total - reads_unique
    
    @property
    def df(self):
        df = pd.DataFrame()
        df['stat'] = [self.reads_total, self.reads_unique, self.reads_removed]
        df['stat_type'] = ['reads_total', 'reads_unique', 'reads_removed']
        df['read_type'] = self.read_type
        df['read_number'] = 1
        df['sample'] = self.sample
        return df

class ReadDeduplicationProcess(Process):
    def __init__(
        self,
        inq: Queue,
        outq: Queue,
        hash_seed: int = 42,
        save_hashed_dict_path: str = "",
    ):

        self.inq = inq
        self.outq = outq
        self.hash_seed = hash_seed
        self.save_hashed_dict_path = save_hashed_dict_path

        super(ReadDeduplicationProcess, self).__init__()

    def _save_dict(self, d):
        if self.save_hashed_dict_path:
            with open(self.save_hashed_dict_path, "wb") as w:
                pickle.dump(d, w)

    def run(self):
        hash_seed = self.hash_seed
        deduplicated_dict = dict()

        reads = self.inq.get()
        while not reads == "END":

            for read_glob in reads:
                hash_sequence = mmh3.hash64(
                    "".join([r.sequence for r in read_glob]), seed=hash_seed
                )[0]
                hash_id = mmh3.hash64(
                    "".join([r.name for r in read_glob]), seed=hash_seed
                )[0]

                if not hash_sequence in deduplicated_dict:
                    deduplicated_dict[hash_sequence] = hash_id

            reads = self.inq.get()

        self._save_dict(deduplicated_dict)
        self.outq.put("END")


class ReadDuplicateRemovalProcess(Process):
    def __init__(
        self,
        inq: Queue,
        outq: Queue,
        deduplicated_ids: set,
        statq: Queue = None,
        hash_seed: int = 42,
    ):

        self.inq = inq
        self.outq = outq
        self.hash_seed = hash_seed
        self.deduplicated_ids = deduplicated_ids
        self.statq = statq
        self.reads_total = 0
        self.reads_unique = 0

        super(ReadDuplicateRemovalProcess, self).__init__()

    def run(self):

        hash_seed = self.hash_seed
        deduplicated_ids = self.deduplicated_ids
        reads_unique = list()

        reads = self.inq.get()
        while not reads == "END":
            for read_glob in reads:
                hash_id = mmh3.hash64(
                    "".join([r.name for r in read_glob]), seed=hash_seed
                )[0]

                if hash_id in deduplicated_ids:
                    print(hash_id)
                    reads_unique.append(read_glob)

            self.outq.put(reads_unique)

            self.reads_total += len(reads)
            self.reads_unique += len(reads_unique)
            reads_unique = list()
            reads = self.inq.get()

        self.outq.put("END")
    
        if self.statq:
            self.statq.put({'reads_total': self.reads_total, 'reads_unique': self.reads_unique})
            self.statq.put('END')


def main(
    mode: str,
    input_files: list,
    output_files: list = None,
    deduplicated_ids: Union[str, list] = None,
    read_buffer=10000,
    compression_level=5,
    n_cores=1,
):

    # Set up multiprocessing variables
    inputq = SimpleQueue()  # Reads are placed into this queue for deduplication
    writeq = SimpleQueue()  # Deduplicated reads are placed into the queue for writing
    statq = SimpleQueue()   # Statistics are sent on this queue for processing

    if mode == "find_duplicates":
        reader = FastqReaderProcess(
            input_files=input_files,
            outq=inputq,
            n_subprocesses=n_cores,
            read_buffer=read_buffer,
        )

        deduplicator = ReadDeduplicationProcess(
            inq=inputq, outq=writeq, save_hashed_dict_path=deduplicated_ids
        )

        processes = [reader, deduplicator]
        
        for proc in processes:
            proc.start()

        for proc in processes:
            proc.join()
            proc.terminate()
    
    elif mode == 'merge_ids':

        dicts = []
        for dd_id in input_files:
            with open(dd_id, "rb") as r:
                dicts.append(pickle.load(r))

        dedup_dict = merge_dictionaries(dicts)  # {SEQUENCE_HASH: READ_NAME_HASH}

        with open(output_files, 'wb') as w:
            pickle.dump(dedup_dict, w)
        

    elif mode == "remove_duplicates":

        dicts = []
        for dd_id in deduplicated_ids:
            with open(dd_id, "rb") as r:
                dicts.append(pickle.load(r))

        dedup_dict = merge_dictionaries(dicts)  # {SEQUENCE_HASH: READ_NAME_HASH}
        deduplicated_ids_set = set(
            dedup_dict.values()
        )  # {READ_NAME_HASH_1, READ_NAME_HASH_2}

        reader = FastqReaderProcess(
            input_files=input_files,
            outq=inputq,
            read_buffer=read_buffer,
            n_subprocesses=n_cores,
        )

        deduplicator = [
            ReadDuplicateRemovalProcess(
                inq=inputq, outq=writeq, deduplicated_ids=deduplicated_ids_set, statq=statq
            )
            for _ in range(n_cores)
        ]

        writer = FastqWriterProcess(
            inq=writeq,
            output=output_files,
            compression_level=compression_level,
        )

        processes = [reader, *deduplicator, writer]
        
        for proc in processes:
            proc.start()

        for proc in processes:
            proc.join()
            proc.terminate()
        
        # Handle statistics
        stats_aggregator = Counter()
        for i in range(n_cores):
            stats_aggregator.update(statq.get())
        
        stats_file_prefix = input_files[0].split('.fastq')[0].rstrip('_12')
        deduplication_stats = DeduplicationStatistics(sample=stats_file_prefix, **stats_aggregator)
        print(deduplication_stats.df)
        deduplication_stats.df.to_csv(f'{stats_file_prefix}.csv')
