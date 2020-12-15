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
from ccanalyser.utils.tools import DigestedRead
from multiprocessing import SimpleQueue, Process

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
                    stats_dict = self._append_statistics(fragment_stats, stat_dict)
          
        

        return self._format_stats_dict(stat_dict)
  
class DigestionStatistics():
    def __init__(self,
                 sample: str,
                 read_type='pe',
                 read_number=1,  
                 slices_total: Counter = None,
                 slices_valid: Counter = None):

        self.sample = sample
        self.read_type = read_type
        self.read_number = read_number
        self.slices_total = slices_total
        self.slices_valid = slices_valid
    
    @property
    def n_total_slices(self):
        return sum(n_slices * count for n_slices, count in self.slices_total.items())
    
    @property
    def n_total_reads(self):
        return sum(self.slices_total.values())
    
    @property
    def n_valid_slices(self):
        return sum(n_slices * count for n_slices, count in self.slices_valid.items())
    
    @property
    def n_valid_reads(self):
        return sum(self.slices_valid.values())
    
    @property
    def unfiltered_histogram(self):
        df = pd.DataFrame()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['number_of_slices'] = self.slices_total.keys()
        df['count'] = self.slices_total.values()
        return df
    
    @property
    def filtered_histogram(self):
        df = pd.DataFrame()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['number_of_slices'] = self.slices_valid.keys()
        df['count'] = self.slices_valid.values()
        return df
    
    @property
    def slice_summary_df(self):
        df = pd.DataFrame()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['stat_type'] = ['total', 'valid']
        df['stat'] = [self.n_total_slices, self.n_valid_slices]
        return df
    
    @property
    def read_summary_df(self):
        df = pd.DataFrame()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['stat_type'] = ['total', 'valid']
        df['stat'] = [self.n_total_reads, self.n_valid_reads]
        return df
    
    def __add__(self, other):
        self.slices_total.update(other.slices_total)
        self.slices_valid.update(other.slices_valid)

        return self


class ReadDigestionProcess(Process):
    def __init__(
        self,
        inq: SimpleQueue,
        outq: SimpleQueue,
        statq: SimpleQueue = None,
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
                ].slices_valid_counter
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

                if all(digested_str):  # Make sure that all reads have valid slices
                    buffer_reads.append("".join(digested_str))

                
                digestion_stats = {read_number + 1: {'total': d.slices_total_counter, 'valid': d.slices_valid}
                                   for read_number, d in enumerate(digested)}
                
                buffer_stats.append(digestion_stats)

            self.outq.put("".join(buffer_reads))

            if self.statq:
                self.statq.put(buffer_stats)

            buffer_reads = []
            buffer_stats = []
            reads = self.inq.get()

        self.outq.put("END")

        if self.statq:
            self.statq.put('END')









def main(
    subcommand,
    input_fastq=None,
    fq1=None,
    fq2=None,
    restriction_enzyme=None,
    keep_cutsite=False,
    output_file="out.fastq.gz",
    minimum_slice_length=18,
    stats_file=None,
    compression_level=5,
    n_digestion_processes=1,
    buffer=10000,
):

    # Set up multiprocessing variables
    inputq = SimpleQueue()  # reads are placed into this queue for processing
    writeq = SimpleQueue()  # digested reads are placed into the queue for writing
    statq = SimpleQueue()  # stats queue

    # Variables
    cut_site = get_re_site(restriction_enzyme)

    if subcommand == "flashed":  # Checks the subsubcommand to see in which mode to run

        fq_reader = FastqReaderProcess(
            input_files=input_fastq, outq=inputq, n_subprocesses=n_digestion_processes
        )

        digestion_processes = [
            ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=subcommand,
            )
            for _ in range(n_digestion_processes)
        ]

    elif subcommand == "unflashed":

        fq_reader = FastqReaderProcess(
            input_files=[fq1, fq2], outq=inputq, n_subprocesses=n_digestion_processes
        )

        digestion_processes = [
            ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=subcommand,
                allow_undigested=True,
            )
            for _ in range(n_digestion_processes)
        ]

    # Writer process is common to both
    fq_writer = FastqWriterProcess(
        inq=writeq, output=output_file, n_subprocesses=n_digestion_processes
    )

    # Start all processes
    processes = [fq_writer, fq_reader, *digestion_processes]

    for proc in processes:
        proc.start()

    # Join processes (wait for processes to finish before the main process)
    fq_writer.join()

    for proc in processes:
        proc.join()
        proc.terminate()
    
    # Collate stats
    collated_stats = DigestionStatCollector(statq, n_digestion_processes).get_collated_stats()

    for read_number, stats in collated_stats.items():
        digestion_stats = DigestionStatistics(output_file,
                                             read_type=subcommand,
                                             read_number=read_number,
                                             slices_total=stats['total'],
                                             slices_valid=stats['valid'])


        
        digestion_stats.unfiltered_histogram.to_csv(f'{stats_file}.unfiltered.histogram.tsv')
        digestion_stats.filtered_histogram.to_csv(f'{stats_file}.filtered.histogram.tsv')
        digestion_stats.slice_summary_df.to_csv(f'{stats_file}.slice.summary.tsv')
        digestion_stats.read_summary_df.to_csv(f'{stats_file}.read.summary.tsv')

        

    



    
    
