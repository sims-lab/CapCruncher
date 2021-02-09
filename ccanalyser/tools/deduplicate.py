import ujson
import sys
from multiprocessing import Process, Queue, SimpleQueue
from xopen import xopen
import xxhash


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

            self.outq.put(reads_unique)

            self.reads_total += len(reads)
            self.reads_unique += len(reads_unique)
            reads_unique = list()
            reads = self.inq.get()

        self.outq.put("END")
    
        if self.statq:
            self.statq.put({'reads_total': self.reads_total, 'reads_unique': self.reads_unique})
            self.statq.put('END')



