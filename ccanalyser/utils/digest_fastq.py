from collections import Counter, defaultdict
from multiprocessing.queues import Queue
from os import stat
from typing import Union
from numpy.lib.histograms import histogram

import pandas as pd
import numpy as np
from pysam import FastxFile
from xopen import xopen

from ccanalyser.utils.helpers import get_re_site
from ccanalyser.utils.io import FastqReaderProcess, FastqWriterProcess
from multiprocessing import SimpleQueue, Process, Queue
import re
from numpy.random import randint

class DigestionStatCollector():
    def __init__(self, statq: Queue, n_subprocess: int = 1) -> None:

        self.statq = statq
        self.n_subprocesses = n_subprocess
        self.n_subprocesses_terminated = 0
        self.stat_dict = None
    
    def _set_up_stats_dict(self, stat_sample):
        
        stat_dict = dict()
        
        for rn in stat_sample:
            stat_dict[rn] = defaultdict(list)
        
        return stat_dict

    def _append_statistics(self, frag_stats, stats_dict):
        
        for read_number, read_stats in frag_stats.items():
            for stat_type, count in read_stats.items():
                stats_dict[read_number][stat_type].append(count)
        
        return stats_dict


    def _format_stats_dict(self, stats_dict):
        for read_number, read_number_stats in stats_dict.items():
            for stat_type, counts in read_number_stats.items():
                stats_dict[read_number][stat_type] = Counter(counts)
        
        return stats_dict


    def get_collated_stats(self) -> dict:
              
        stats = self.statq.get()
        stat_dict = self._set_up_stats_dict(stats[0])
        
        while self.n_subprocesses_terminated < self.n_subprocesses:
            
           
            if stats == 'END':
                self.n_subprocesses_terminated += 1
                continue
        
            else:
                for fragment_stats in stats:
                    stat_dict = self._append_statistics(fragment_stats, stat_dict)
            
            stats = self.statq.get()

        return self._format_stats_dict(stat_dict)
  
class DigestionStatistics():
    def __init__(self,
                 sample: str,
                 read_type='pe',
                 read_number=1,  
                 slices_unfiltered: Counter = None,
                 slices_filtered: Counter = None):

        self.sample = sample
        self.read_type = read_type
        self.read_number = read_number
        self.slices_unfiltered = slices_unfiltered
        self.slices_filtered = slices_filtered


        self.n_unfiltered_slices = self._get_n_unfiltered_slices()
        self.n_unfiltered_reads = self._get_n_unfiltered_reads()
        self.n_filtered_slices = self._get_n_filtered_slices()
        self.n_filtered_reads = self._get_n_filtered_reads()
        self.filtered_histogram = self._get_filtered_histogram()
        self.unfiltered_histogram = self._get_unfiltered_histogram()
        self.slice_summary = self._get_slice_summary_df()
        self.read_summary = self._get_read_summary_df()

  

    def _get_n_unfiltered_slices(self):
        return sum(n_slices * count for n_slices, count in self.slices_unfiltered.items())
    
    def _get_n_unfiltered_reads(self):
        return sum(self.slices_unfiltered.values())
    
    def _get_n_filtered_slices(self):
        return sum(n_slices * count for n_slices, count in self.slices_filtered.items())

    def _get_n_filtered_reads(self):
        return sum(v for k, v in self.slices_filtered.items() if not k == 0)
    
    def _get_unfiltered_histogram(self):
        df = pd.DataFrame()
        df['number_of_slices'] = self.slices_unfiltered.keys()
        df['count'] = self.slices_unfiltered.values()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        return df[['sample', 'read_type', 'read_number', 'number_of_slices', 'count']].sort_values('number_of_slices')
    
    def _get_filtered_histogram(self):
        df = pd.DataFrame()
        df['number_of_slices'] = self.slices_filtered.keys()
        df['count'] = self.slices_filtered.values()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        return df[['sample', 'read_type', 'read_number', 'number_of_slices', 'count']].sort_values('number_of_slices')
    
    def _get_slice_summary_df(self):
        df = pd.DataFrame()
        df['stat_type'] = ['unfiltered', 'filtered']
        df['stat'] = [self.n_unfiltered_slices, self.n_filtered_slices]
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['stage'] = 'digestion'
        return df[['sample','stage', 'read_type', 'read_number', 'stat_type', 'stat']]
    
    def _get_read_summary_df(self):
        df = pd.DataFrame()
        df['stat_type'] = ['unfiltered', 'filtered']
        df['stat'] = [self.n_unfiltered_reads, self.n_filtered_reads]
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['stage'] = 'digestion'
        return df[['sample', 'stage','read_type', 'read_number', 'stat_type', 'stat']]
    
    def __add__(self, other):       
        self.unfiltered_histogram = pd.concat([self.unfiltered_histogram, other.unfiltered_histogram]).reset_index(drop=True)
        self.filtered_histogram = pd.concat([self.filtered_histogram, other.filtered_histogram]).reset_index(drop=True)
        self.slice_summary = pd.concat([self.slice_summary, other.slice_summary]).reset_index(drop=True)
        self.read_summary = pd.concat([self.read_summary, other.read_summary]).reset_index(drop=True)
        return self

class DigestedRead:
    def __init__(
        self,
        read,
        cutsite,
        min_slice_length=18,
        slice_number_offset=0,
        allow_undigested=False,
        read_type="flashed",
    ):

        self.read = read
        self.min_slice_length = min_slice_length
        self.slice_number_offset = slice_number_offset
        self.allow_undigested = allow_undigested
        self.read_type = read_type

        self.recognition_seq = cutsite.upper()
        self.recognition_len = len(cutsite)
        self.recognition_re = re.compile(self.recognition_seq)

        self.slice_indexes = self.get_recognition_site_indexes()
        self.slices_unfiltered = len(self.slice_indexes) - 1
        self.slices_filtered = 0
        self.has_slices = self.slices_unfiltered > 1 
        self.slices = self._get_slices()

    def get_recognition_site_indexes(self):
        indexes = [
            re_site.start()
            for re_site in self.recognition_re.finditer(self.read.sequence.upper())
        ]

        indexes.insert(0, 0)
        indexes.append(len(self.read.sequence))

        return indexes

    def _get_slices(self):

        indexes = self.slice_indexes
        slice_no = self.slice_number_offset
        slices_list = []

        if self.has_slices or self.allow_undigested:

            # Iterate through offset indexes to get correct start and end
            for ii, (slice_start, slice_end) in enumerate(zip(indexes, indexes[1:])):

                # If this is not the first slice
                if ii > 0:
                    slice_start += self.recognition_len

                if self.is_filtered_slice(slice_start, slice_end):
                    slices_list.append(
                        self.prepare_slice(slice_start, slice_end, slice_no)
                    )

                    self.slices_filtered += 1
                    slice_no += 1

        return slices_list

    def prepare_slice(self, start, end, slice_no):
        return "\n".join(
            [
                f"@{self.read.name}|{self.read_type}|{slice_no}|{randint(0,100)}",
                self.read.sequence[start:end],
                "+",
                self.read.quality[start:end],
            ]
        )

    def is_filtered_slice(self, start, end):
        if (end - start) >= self.min_slice_length:
            return True
    
    def __str__(self):
        return ("\n".join(self.slices) + "\n") if self.slices else ""

class ReadDigestionProcess(Process):
    def __init__(
        self,
        inq: SimpleQueue,
        outq: SimpleQueue,
        statq: Queue = None,
        **digestion_kwargs
    ) -> None:

        super(ReadDigestionProcess, self).__init__()

        self.inq = inq
        self.outq = outq
        self.statq = statq
        self.digestion_kwargs = digestion_kwargs
        self.read_type = digestion_kwargs["read_type"]

    def _digest_reads(self, reads, **digestion_kwargs):
        digested = []
        for i, read in enumerate(reads):
            if i == 0:
                digested.append(DigestedRead(read, **digestion_kwargs))
            else:
                digestion_kwargs["slice_number_offset"] = digested[
                    i - 1
                ].slices_filtered
                digested.append(DigestedRead(read, **digestion_kwargs))

        return digested

    def run(self):

        reads = self.inq.get()
        buffer_reads = []
        buffer_stats = []

        while not reads == "END":

            for read in reads:
                digested = self._digest_reads(read, **self.digestion_kwargs)
                digested_str = [str(dr) for dr in digested]
                digestion_stats = {read_number + 1: {'unfiltered': d.slices_unfiltered, 'filtered': d.slices_filtered}
                                   for read_number, d in enumerate(digested)}
                buffer_stats.append(digestion_stats)

                if all(digested_str):  # Make sure that all reads have filtered slices
                    buffer_reads.append("".join(digested_str))

            self.outq.put("".join(buffer_reads))

            if self.statq:
                self.statq.put_nowait(buffer_stats)

            buffer_reads = []
            buffer_stats = []
            reads = self.inq.get()

        self.outq.put("END")

        if self.statq:
            self.statq.put_nowait('END')
        
def main(
    subcommand,
    input_fastq=None,
    restriction_enzyme=None,
    keep_cutsite=False,
    output_file="out.fastq.gz",
    minimum_slice_length=18,
    compression_level=5,
    n_digestion_processes=1,
    buffer=100000,
    stats_prefix='',
):

    # Set up multiprocessing variables
    inputq = SimpleQueue()  # reads are placed into this queue for processing
    writeq = SimpleQueue()  # digested reads are placed into the queue for writing
    statq = Queue()  # stats queue

    # Variables
    cut_site = get_re_site(restriction_enzyme)

    if subcommand == "flashed":  # Checks the subsubcommand to see in which mode to run

        reader = FastqReaderProcess(
            input_files=input_fastq, outq=inputq, n_subprocesses=n_digestion_processes, read_buffer=buffer,
        )

        digestion_processes = [
            ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=subcommand,
                allow_undigested=False,
                statq=statq,
            )
            for _ in range(n_digestion_processes)
        ]

    elif subcommand == "pe":

        reader = FastqReaderProcess(
            input_files=input_fastq, outq=inputq, n_subprocesses=n_digestion_processes, read_buffer=buffer
        )

        digestion_processes = [
            ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=subcommand,
                allow_undigested=True,
                statq=statq,
            )
            for _ in range(n_digestion_processes)
        ]

    # Writer process is common to both
    writer = FastqWriterProcess(
        inq=writeq, output=output_file, n_subprocesses=n_digestion_processes, compression_level=compression_level
    )

    # Start all processes
    processes = [writer, reader, *digestion_processes]

    for proc in processes:
        proc.start()


    reader.join()
    writer.join()
    
    # Collate stats
    print('')
    print('Collating stats')
    collated_stats = DigestionStatCollector(statq, n_digestion_processes).get_collated_stats()
    sample_name = re.match(r'.*/(.*)(_part\d+.).*', input_fastq[0]).group(1)
    
    stats = [DigestionStatistics(sample=sample_name,
                                 read_type=subcommand,
                                 read_number=read_number,
                                 slices_unfiltered=stats['unfiltered'],
                                 slices_filtered=stats['filtered'])
                        for read_number, stats in collated_stats.items()
                        ]
    
    if len(stats) > 1: # Need to collate stats from digestion of 2+ files
        for stat in stats[1:]:
            digestion_stats = stats[0] + stat
    else:
        digestion_stats = stats[0]


    digestion_stats.unfiltered_histogram.to_csv(f'{stats_prefix}.digestion.unfiltered.histogram.csv')
    digestion_stats.filtered_histogram.to_csv(f'{stats_prefix}.digestion.filtered.histogram.csv')
    digestion_stats.slice_summary.to_csv(f'{stats_prefix}.digestion.slice.summary.csv')
    digestion_stats.read_summary.to_csv(f'{stats_prefix}.digestion.read.summary.csv')

    print(digestion_stats.read_summary)


    



    
    
