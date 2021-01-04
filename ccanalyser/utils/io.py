import pathlib
from multiprocessing import Manager, Process, Queue
from typing import Union

import numpy as np
import pandas as pd
from pysam import FastxFile
from xopen import xopen


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

        while self.n_workers_terminated < self.n_subprocesses:
            

            if reads == 'END':
                self.n_workers_terminated += 1
                continue
            
            elif is_string_input:
                self.file_handles[0].write(reads)
            
            else:
                reads_str = ['\n'.join([str(r) for r in read_glob]) for read_glob in zip(*reads)]
                
                for fh, read_set in zip(self.file_handles, reads_str):
                    fh.write((read_set + '\n'))
            
            reads = self.inq.get()
        
        for fh in self.file_handles:
            fh.close()
                




class FastqWriterSplitterProcess(Process):
    def __init__(self,
                 inq: Queue,
                 output_prefix: Union[str, list],
                 compression_level: int = 3,
                 compression_threads: int = 8,
                 n_subprocesses: int = 1,
                 paired_output=False,
                 ):
        

        self.inq = inq
        self.output_prefix = output_prefix
        self.paired_output = paired_output

        self.compression_level = compression_level
        self.compression_threads = compression_threads

        self.n_subprocesses = n_subprocesses
        self.n_workers_terminated = 0
        
        super(FastqWriterSplitterProcess, self).__init__()

    def run(self):

        n_chunk = 0
        while self.n_workers_terminated < self.n_subprocesses:
           
            reads = self.inq.get()
           
            if reads == 'END':
                self.n_workers_terminated += 1
                continue
            else:
                if self.paired_output:
                    r1, r2 = zip(*((str(r1), str(r2)) for r1, r2 in reads))
                    fn1 = f'{self.output_prefix}_part{n_chunk}_1.fastq.gz'
                    fn2 = f'{self.output_prefix}_part{n_chunk}_2.fastq.gz'

                    with xopen(fn1, 'w', compresslevel=self.compression_level, threads=self.compression_threads // 2) as w1:
                        with xopen(fn2, 'w', compresslevel=self.compression_level, threads=self.compression_threads // 2) as w2:
                            w1.write(('\n'.join(r1) + '\n'))
                            w2.write(('\n'.join(r2) + '\n'))
                    
                else:
                    r = [str(r) for r in reads]
                    fn = f'{self.output_prefix}_part{n_chunk}.fastq.gz'
                    with xopen(fn, 'w', compresslevel=self.compression_level, threads=1) as w:
                        w.write(('\n'.join(r) + '\n'))
        
            n_chunk += 1
            







# class CollateDigestionStatsProcess(Process):
#     def __init__(self,
#                  inq: Queue,
#                  n_subprocesses: int = 1,
#                  digestion_mode=None) -> None:
        
#         super(ReadDigestionProcess, self).__init__(self)

#         self.inq = inq
#         self.n_workers_terminated = 0
#         self.n_subprocesses = n_subprocesses
#         self.digestion_mode = digestion_mode

#     def run(self):

#         stats = []
#         counts = self.inq.get()



#         while self.n_workers_terminated < self.n_subprocesses:

     
#             if counts == "END":
#                 self.n_workers_terminated += 1
#                 continue
            
                
            
        

        


        

