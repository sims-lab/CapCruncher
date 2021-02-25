import re
import pysam
from multiprocessing import Queue, Process, SimpleQueue
from numpy.random import randint
from typing import Union
import traceback

class DigestedChrom:
    def __init__(
        self,
        chrom: pysam.FastqProxy,
        cutsite: str,
        fragment_number_offset: int = 0,
        fragment_min_len: int = 1,
    ):

        self.chrom = chrom
        self.recognition_seq = cutsite.upper()
        self.recognition_len = len(cutsite)
        self.recognition_re = re.compile(self.recognition_seq)

        self.fragment_indexes = self.get_recognition_site_indexes()
        self.fragment_number_offset = fragment_number_offset
        self.fragment_min_len = fragment_min_len

    def get_recognition_site_indexes(self):
        indexes = [
            re_site.start()
            for re_site in self.recognition_re.finditer(self.chrom.sequence.upper())
        ]

        indexes.insert(0, 0)
        indexes.append(len(self.chrom.sequence))

        return indexes
    
    @property
    def fragments(self):

        indexes = self.fragment_indexes
        fragment_no = self.fragment_number_offset

        # Iterate through offset indexes to get correct start and end
        for ii, (fragment_start, fragment_end) in enumerate(zip(indexes, indexes[1:])):

            # If this is not the first fragment
            if ii > 0:
                fragment_start += self.recognition_len
            
            # Check to see if the fragment is long enough to be recorded (default 1bp)
            if (fragment_end - fragment_start) >= self.fragment_min_len:
                yield self._prepare_fragment(fragment_start, fragment_end, fragment_no)
                fragment_no += 1

    def _prepare_fragment(self, start, end, fragment_no):
        return '\t'.join([str(x) for x in (self.chrom.name, start, end, fragment_no)]) + '\n'

class DigestedRead:
    def __init__(
        self,
        read, #pysam proxy object
        cutsite: str,
        min_slice_length: int = 18,
        slice_number_offset: int = 0,
        allow_undigested: bool = False,
        read_type: str = "flashed",
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
        self.read_type = digestion_kwargs.get("read_type", 'flashed')

        
        if not 'cutsite' in digestion_kwargs:
            raise KeyError('Cutsite is required to be present in digestion arguments')
            traceback.format_exc()
            self.outq.put('END')


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
        try:
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
        
        except Exception as e:
            traceback.format_exc()
            self.outq.put('END')

  