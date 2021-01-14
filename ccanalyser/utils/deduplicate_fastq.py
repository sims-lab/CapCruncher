#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""
import ujson
import sys
from multiprocessing import Process, Queue, SimpleQueue
from typing import Union
import pandas as pd
from collections import Counter, namedtuple
import xxhash
import re
from ccanalyser.utils.io import FastqReaderProcess, FastqWriterProcess
from xopen import xopen
import os

def open_logfile(fn):
    if not isinstance(fn, type(sys.stdout)):
        return xopen(fn, "w")
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
        df['read_number'] = 0
        df['stage'] = 'deduplication'
        df['sample'] = self.sample
        return df

class ReadDeduplicationParserProcess(Process):
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
        self.read_data = dict()

        super(ReadDeduplicationParserProcess, self).__init__()

    def _save_dict(self, d):
        if self.save_hashed_dict_path:
            with xopen(self.save_hashed_dict_path, "w") as w:
                ujson.dump(d, w)

    def run(self):

        hash_seed = self.hash_seed
        hash_function = xxhash.xxh64_intdigest
        reads = self.inq.get()
        read_data = self.read_data

        while not reads == "END":

            for read_glob in reads:
                hash_sequence = hash_function(
                    "".join([r.sequence for r in read_glob]))
                hash_id = hash_function(
                    "".join([r.name for r in read_glob]))

                read_data[hash_id] = hash_sequence
                #read_data["".join([r.name for r in read_glob])] = "".join([r.sequence for r in read_glob])

            reads = self.inq.get()

        self._save_dict(read_data)
        self.outq.put("END")


class ReadDuplicateRemovalProcess(Process):
    def __init__(
        self,
        inq: Queue,
        outq: Queue,
        duplicated_ids: set,
        statq: Queue = None,
        hash_seed: int = 42,
    ):

        self.inq = inq
        self.outq = outq
        self.hash_seed = hash_seed
        self.duplicated_ids = duplicated_ids
        self.statq = statq
        self.reads_total = 0
        self.reads_unique = 0

        super(ReadDuplicateRemovalProcess, self).__init__()

    def run(self):

        hash_seed = self.hash_seed
        hash_function = xxhash.xxh64_intdigest
        duplicated_ids = self.duplicated_ids
        reads_unique = list()

        reads = self.inq.get()
        while not reads == "END":
            for read_glob in reads:

                hash_id = hash_function(
                     "".join([r.name for r in read_glob]))

                #hash_id = "".join([r.name for r in read_glob])

                if not str(hash_id) in duplicated_ids:
                    reads_unique.append(read_glob)
                else:
                    print(f'{hash_id} matches ')

            self.outq.put(reads_unique)

            self.reads_total += len(reads)
            self.reads_unique += len(reads_unique)
            reads_unique = list()
            reads = self.inq.get()

        self.outq.put("END")
    
        if self.statq:
            self.statq.put({'reads_total': self.reads_total, 'reads_unique': self.reads_unique})
            self.statq.put('END')


def load_json(fn):
    with xopen(fn) as r:
        d = ujson.load(r)
        return d


def main(
    mode: str,
    input_files: list,
    output_files: list = None,
    read_ids: Union[str, list] = None,
    read_buffer=100000,
    compression_level=5,
    n_cores=1,
    stats_prefix: str = None,
    sample_name: str = None,
):

    # Set up multiprocessing variables
    inputq = SimpleQueue()  # Reads are placed into this queue for deduplication
    writeq = SimpleQueue()  # Deduplicated reads are placed into the queue for writing
    statq = SimpleQueue()   # Statistics are sent on this queue for processing

    if mode == "parse":
        reader = FastqReaderProcess(
            input_files=input_files,
            outq=inputq,
            n_subprocesses=n_cores,
            read_buffer=read_buffer,
        )

        parser = ReadDeduplicationParserProcess(
            inq=inputq, outq=writeq, save_hashed_dict_path=read_ids
        )

        processes = [reader, parser]
        
        for proc in processes:
            proc.start()

        for proc in processes:
            proc.join()
            proc.terminate()
    
    elif mode == 'identify':
        
        dedup_sequences = dict()
        read_ids = set()
        
        for fn in input_files:
            d = load_json(fn)
            read_ids.update(d)
            dedup_sequences.update(invert_dictionary(d)) # {READ_NAME_HASH: SEQUENCE_HASH} -> {SEQUENCE_HASH: READ_NAME_HASH}
        
        duplicated_ids = read_ids - set(dedup_sequences.values())

        with xopen(output_files, 'w') as w:
            duplicated_ids_dict = dict.fromkeys(duplicated_ids)
            ujson.dump(duplicated_ids_dict, w)

    elif mode == "remove":

        duplicated_ids = set(load_json(read_ids))

        print(duplicated_ids)

        reader = FastqReaderProcess(
            input_files=input_files,
            outq=inputq,
            read_buffer=read_buffer,
            n_subprocesses=n_cores,
        )

        writer = FastqWriterProcess(
            inq=writeq,
            output=output_files,
            compression_level=compression_level,
        )

        deduplicator = [
            ReadDuplicateRemovalProcess(
                inq=inputq, outq=writeq, duplicated_ids=duplicated_ids, statq=statq
            )
            for _ in range(n_cores)
        ]

        reader.start()
        writer.start()
        for dedup in deduplicator:
            dedup.start()

        processes = [writer, reader, *deduplicator]

        for proc in processes:
            proc.join()
            proc.terminate()
        
        # Handle statistics
        stats_aggregator = Counter()
        stats = statq.get()

        while not stats == 'END':
            stats_aggregator.update(stats)
            stats = statq.get() 


        deduplication_stats = DeduplicationStatistics(sample=sample_name, **stats_aggregator)
        print(deduplication_stats.df)
        deduplication_stats.df.to_csv(f'{stats_prefix}.deduplication.csv')
