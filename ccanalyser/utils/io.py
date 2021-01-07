from multiprocessing.queues import SimpleQueue
import pathlib
from multiprocessing import Manager, Process, Queue, Pipe
from typing import Union
import traceback

import numpy as np
import pandas as pd
from pysam import FastxFile
from xopen import xopen

#TODO: Implement error queues/pipes for handling exceptions

class FastqReaderProcess(Process):
    def __init__(self,
                 input_files: Union[str, list],
                 outq: Queue,
                 read_buffer: int = 100000,
                 read_counter: Manager().Value = None,
                 n_subprocesses: int = 1,
                 statq: Queue = None) -> None:
        
        
        # Input variables
        self.input_files = input_files
        self.multifile = self._is_multifile(input_files)

        if self.multifile:
            self.input_files_pysam = [FastxFile(f) for f in self.input_files]
        else:
            self.input_files_pysam  = [FastxFile(self.input_files), ]
        
        # Multiprocessing variables
        self.outq = outq
        self.statq = statq
        self.n_subprocesses = n_subprocesses

       
        # Reader variables
        self.read_buffer = read_buffer
        self.read_counter = read_counter

        super(FastqReaderProcess, self).__init__()
    
    def _is_multifile(self, files):
        if not isinstance(files, (str, pathlib.Path)):
            return True 
        elif isinstance(files, (list, tuple)) and len(files > 1):
            return True 
        else:
            return False


    def run(self):

        try:
            buffer = []
            for read_counter, read in enumerate(zip(*self.input_files_pysam)):

                buffer.append(read)

                if read_counter % self.read_buffer == 0 and not read_counter == 0:
                    self.outq.put(buffer)
                    buffer = []
                    print(f'Read {read_counter} reads')
            
            self.outq.put(buffer)
            print(f'Number of reads processed: {read_counter + 1}')

            for i in range(self.n_subprocesses):
                self.outq.put('END')

            # Deal with number of reads that have been read
            if self.read_counter:
                self.read_counter.value = read_counter
            elif self.statq:
                self.statq.put({'reads_total': read_counter})
        
        except Exception as e:
            print(traceback.format_exc())
            self.outq.put('END')


class FastqReadFormatterProcess(Process):
    def __init__(self,
                 inq: SimpleQueue,
                 outq: SimpleQueue,
                 formatting: list = None) -> None:
        
        self.inq = inq
        self.outq = outq
        self.formatting = [self._format_as_str, ] if not formatting else formatting

        super(FastqReadFormatterProcess, self).__init__()

    def _format_as_str(self, reads):

        # [(r1, r2), (r1, r2)] -> [r1 combined string, r2 combined string]
        return ['\n'.join([str(rn) for rn in r]) for r in zip(*reads)]

    
    def run(self):

        reads = self.inq.get()
        
        while not reads == "END":
            for formatting_to_apply in self.formatting:
                reads = formatting_to_apply(reads)

            self.outq.put(reads)
            reads = self.inq.get()
        
        
        self.outq.put('END')


class FastqWriterSplitterProcess(Process):
    def __init__(self,
                 inq: Queue,
                 output_prefix: Union[str, list],
                 paired_output: bool = False,
                 compression_level: int = 3,
                 compression_threads: int = 8,
                 n_subprocesses: int = 1,
                 n_workers_terminated: int = 0,
                 n_files_written: int = 0,
                 ):
        

        self.inq = inq
        self.output_prefix = output_prefix
        self.paired_output = paired_output

        self.compression_level = compression_level
        self.compression_threads = compression_threads

        self.n_subprocesses = n_subprocesses
        self.n_workers_terminated = n_workers_terminated
        self.n_files_written = n_files_written
        
        super(FastqWriterSplitterProcess, self).__init__()
    
    def _get_file_handles(self):

        if not self.paired_output:
            fnames = [f'{self.output_prefix}_part{self.n_files_written}.fastq.gz', ]
        else:
            fnames = [f'{self.output_prefix}_part{self.n_files_written}_{i+1}.fastq.gz'
                      for i in range(2)]

        return [xopen(fn, 'w', compresslevel=self.compression_level, threads=2) 
                for fn in fnames]
    
    def run(self):


        reads = self.inq.get()
        is_string_input = True if isinstance(reads[0], str) else False

        while self.n_workers_terminated < self.n_subprocesses:

            if reads == 'END':
                self.n_workers_terminated += 1
                continue
            
            elif is_string_input:
                for fh, read in zip(self._get_file_handles(), reads):
                    fh.write(read)
                    fh.close()
            
            else:
                reads_str = ['\n'.join([str(r) for r in read_glob]) for read_glob in zip(*reads)]
                
                for fh, read_set in zip(self._get_file_handles(), reads_str):
                    fh.write((read_set + '\n'))
                    fh.close()
            
            reads = self.inq.get()
            self.n_files_written += 1

class FastqWriterProcess(Process):
    def __init__(self,
                 inq: Queue,
                 output: Union[str, list],
                 compression_level: int = 5,
                 n_subprocesses: int = 1,
                 ):
        
        super(FastqWriterProcess, self).__init__()

        self.inq = inq
        self.output = output
        self.compression_level = compression_level
        self.n_workers_terminated = 0
        self.n_subprocesses = n_subprocesses
        self.file_handles = self._get_filehandles()
        self.name = 'FastqWriter'

    def _get_filehandles(self):
        if isinstance(self.output, str):
            return [xopen(self.output,'w', compresslevel=self.compression_level, threads=2),]
        elif isinstance(self.output, (list, tuple, pd.Series)):
            return [xopen(fn, 'w', compresslevel=self.compression_level, threads=2) for fn in self.output]


    def _inputs_match_number_of_handles(self, reads):
        if len(reads[0]) == len(self.output):
            return True

    def run(self):


        reads = self.inq.get()
        is_string_input = True if isinstance(reads, str) else False

        counter = 0
        while self.n_workers_terminated < self.n_subprocesses:

            if reads == 'END':
                self.n_workers_terminated += 1
                continue
            
            elif is_string_input:
                for fh in self.file_handles:
                    fh.write(reads)
            
            else:
                reads_str = ['\n'.join([str(r) for r in read_glob]) for read_glob in zip(*reads)]
                
                for fh, read_set in zip(self.file_handles, reads_str):
                    fh.write((read_set + '\n'))

            reads = self.inq.get()

        
        for fh in self.file_handles:
            fh.close()
           



             




           
                
            
        

        


        

