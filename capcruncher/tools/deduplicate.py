import queue
from typing import NamedTuple
import ujson
import multiprocessing
from multiprocessing import Process, Queue, SimpleQueue
from xopen import xopen
import xxhash
import functools
import os
import pickle
from collections import namedtuple

#TODO: Look at https://github.com/realead/cykhash/blob/master/doc/README_API.md
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

    def _save_dict(self, d):

        if ".json" in self.output_path:
            with xopen(self.output_path, "w") as w:
                ujson.dump(d, w)
        elif any(ext in self.output_path for ext in [".pkl",".pickle"]):
            with open(self.output_path, "wb") as w:
                pickle.dump(d, w)


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
                        hash_sequence = hash_function("".join([r.sequence for r in read_set]))
                        hash_id = hash_function("".join([r.name for r in read_set]))
                        records[hash_id] = hash_sequence
                
                else:
                    break
            
            except queue.Empty:
                continue

            
        self._save_dict(records)



RemovalStatistics = namedtuple('RemovalStatistics', ["reads_total", "reads_unique", "reads_removed"])


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

        #Misc
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
        
        stats = RemovalStatistics(self.reads_total, self.reads_unique, self.reads_total - self.reads_unique)
        self.stats_tx.send(stats)

