import ujson
import multiprocessing
from multiprocessing import Process, Queue, SimpleQueue
from xopen import xopen
import xxhash
import functools
import os

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
        inq: multiprocessing.SimpleQueue,
        outq: multiprocessing.SimpleQueue,
        hash_seed: int = 42,
        save_hashed_dict_path: os.PathLike = "parsed.json",
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
        """Processes fastq reads from multiple files and generates a hashed json dictionary.
           
           Dictionary is hashed and in the format {(read  1 name + read 2 name): (sequence 1 + sequence 2)}.

           Output path is specified by save_hashed_dict_path.

        """        

        hash_seed = self.hash_seed
        hash_function = functools.partial(xxhash.xxh64_intdigest, seed=hash_seed)
        reads = self.inq.get()
        read_data = self.read_data

        while not reads == "END":

            for read_glob in reads:
                hash_sequence = hash_function("".join([r.sequence for r in read_glob]))
                hash_id = hash_function("".join([r.name for r in read_glob]))

                read_data[hash_id] = hash_sequence
                # read_data["".join([r.name for r in read_glob])] = "".join([r.sequence for r in read_glob])

            reads = self.inq.get()

        self._save_dict(read_data)
        self.outq.put("END")


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
        inq: multiprocessing.SimpleQueue,
        outq: multiprocessing.SimpleQueue,
        duplicated_ids: set,
        statq: multiprocessing.Queue = None,
        hash_seed: int = 42,
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
        self.statq = statq
        self.reads_total = 0
        self.reads_unique = 0

        super(ReadDuplicateRemovalProcess, self).__init__()

    def run(self):

        """Performs read deduplication based on sequence. 

           Unique reads are placed on outq and deduplication statistics are placed on statq.

        """         

        hash_seed = self.hash_seed
        hash_function = functools.partial(xxhash.xxh64_intdigest, seed=hash_seed)
        duplicated_ids = self.duplicated_ids
        reads_unique = list()

        reads = self.inq.get()
        while not reads == "END":
            for read_glob in reads:

                hash_id = hash_function("".join([r.name for r in read_glob]))

                if not hash_id in duplicated_ids:
                    reads_unique.append(read_glob)

            self.outq.put(reads_unique)

            self.reads_total += len(reads)
            self.reads_unique += len(reads_unique)
            reads_unique = list()
            reads = self.inq.get()

        self.outq.put("END")

        if self.statq:
            self.statq.put(
                {"reads_total": self.reads_total, "reads_unique": self.reads_unique}
            )
            self.statq.put("END")
