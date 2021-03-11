import re
import pysam
import multiprocessing
from multiprocessing import Queue, Process, SimpleQueue
import numpy as np
from typing import Iterable, Union, List
import traceback


class DigestedChrom:
    """
    Performs in slico digestion of fasta files.

    Identifies all restriction sites for a supplied restriction enzyme/restriction site
    and generates bed style entries.

    Attributes:
     chrom (pysam.FastqProxy): Chromosome to digest
     recognition_seq (str): Sequence of restriction recognition site
     recognition_len (int): Length of restriction recognition site
     recognition_seq (re.Pattern): Regular expression for restriction recognition site
     fragment_indexes (List[int]): Indexes of fragment(s) start and end positions.
     fragment_number_offset (int): Starting fragment number.
     fragment_min_len (int): Minimum fragment length required to report fragment
           
    """

    def __init__(
        self,
        chrom: pysam.FastqProxy,
        cutsite: str,
        fragment_number_offset: int = 0,
        fragment_min_len: int = 1,
    ):
        """
        Args:
         chrom (pysam.FastqProxy): Input fastq entry to digest
         cutsite (str): Restriction enzyme recognition sequence.
         fragment_number_offset (int, optional): Changes the fragment number output. Useful for multiple digests. Defaults to 0.
         fragment_min_len (int, optional): Minimum length of a fragment required to output. Defaults to 1.
        """        

        self.chrom = chrom
        self.recognition_seq = cutsite.upper()
        self.recognition_len = len(cutsite)
        self.recognition_re = re.compile(self.recognition_seq)

        self.fragment_indexes = self.get_recognition_site_indexes()
        self.fragment_number_offset = fragment_number_offset
        self.fragment_min_len = fragment_min_len

    def get_recognition_site_indexes(self) -> List[int]:
        """
        Gets the start position of all recognition sites present in the sequence.

        Notes:
         Also appends the start and end of the sequnece to enable clearer itteration
         through the indexes.


        Returns:
            List[int]: Indexes of fragment(s) start and end positions.
        """

        indexes = [
            re_site.start()
            for re_site in self.recognition_re.finditer(self.chrom.sequence.upper())
        ]

        indexes.insert(0, 0)
        indexes.append(len(self.chrom.sequence))

        return indexes
    
    @property
    def fragments(self) -> Iterable[str]:
        """
        Extracts the coordinates of restriction fragments from the sequence.

        Yields:
            Iterator[Iterable[str]]: Identified restriction fragments in bed format.
        """

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

    def _prepare_fragment(self, start: int, end: int, fragment_no: int) -> str:
        """Formats fragment into bed style coordinates.

        Args:
         start (int): Fragment start.
         end (int): Fragment end.
         fragment_no (int): Fragment number.

        Returns:
         str: Bed formatted coordinates.
        """
        return '\t'.join([str(x) for x in (self.chrom.name, start, end, fragment_no)]) + '\n'

class DigestedRead:
    """
    Performs in slico digestion of fastq files.

    Identifies all restriction sites for a supplied restriction enzyme/restriction site
    and generates bed style entries.

    Attributes:
     read (pysam.FastqProxy): Read to digest.
     recognition_seq (str): Sequence of restriction recognition site.
     recognition_len (int): Length of restriction recognition site.
     recognition_seq (re.Pattern): Regular expression for restriction recognition site.
     slices (List[str]): List of Fastq formatted digested reads (slices).  
     slice_indexes (List[int]): Indexes of fragment(s) start and end positions.
     slice_number_offset (int): Starting fragment number.
     min_slice_len (int): Minimum fragment length required to report fragment.
     has_slices (bool): Recognition site(s) present within sequence.


           
    """
    def __init__(
        self,
        read: pysam.FastqProxy,
        cutsite: str,
        min_slice_length: int = 18,
        slice_number_offset: int = 0,
        allow_undigested: bool = False,
        read_type: str = "flashed",
    ):
        """
        Args:
         read (pysam.FastqProxy): Read to digest.
         cutsite (str): Restriction enzyme recognition sequence.
         min_slice_length (int, optional): Minimum slice length required for output. Defaults to 18.
         slice_number_offset (int, optional): Increases the reported output slice number. Defaults to 0.
         allow_undigested (bool, optional): If True slices without a restriction site are not filtered out. Defaults to False.
         read_type (str, optional): Combined (flashed) or non-combined (pe). Choose from (flashed|pe). Defaults to "flashed".
        """    

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

    def get_recognition_site_indexes(self) -> List[int]:
        indexes = [
            re_site.start()
            for re_site in self.recognition_re.finditer(self.read.sequence.upper())
        ]

        indexes.insert(0, 0)
        indexes.append(len(self.read.sequence))

        return indexes

    def _get_slices(self) -> List[str]:

        indexes = self.slice_indexes
        slice_no = self.slice_number_offset
        slices_list = []

        if self.has_slices or self.allow_undigested:

            # Iterate through offset indexes to get correct start and end
            for ii, (slice_start, slice_end) in enumerate(zip(indexes, indexes[1:])):

                # If this is not the first slice
                if ii > 0:
                    slice_start += self.recognition_len

                if self._slice_passes_filter(slice_start, slice_end):
                    slices_list.append(
                        self._prepare_slice(slice_start, slice_end, slice_no)
                    )

                    self.slices_filtered += 1
                    slice_no += 1

        return slices_list

    def _prepare_slice(self, start, end, slice_no):
        return "\n".join(
            [
                f"@{self.read.name}|{self.read_type}|{slice_no}|{np.random.randint(0,100)}",
                self.read.sequence[start:end],
                "+",
                self.read.quality[start:end],
            ]
        )

    def _slice_passes_filter(self, start: int, end: int) -> bool:
        """
        Determines if slice exceeds the minimum slice length.

        Args:
         start (int): Slice start position.
         end (int): Slice end position.

        Returns:
         bool: True if greater than minimum slice length  
        
        """    

        if (end - start) >= self.min_slice_length:
            return True
    
    def __str__(self):
        return ("\n".join(self.slices) + "\n") if self.slices else ""

class ReadDigestionProcess(Process):
    """
    Process subclass for multiprocessing fastq digestion.

    """
    def __init__(
        self,
        inq: multiprocessing.SimpleQueue,
        outq: multiprocessing.SimpleQueue,
        statq: multiprocessing.Queue = None,
        **digestion_kwargs
    ) -> None:
        
        """
        Args:
         inq (multiprocessing.SimpleQueue): Queue to hold list of reads to digest.
         outq (multiprocessing.SimpleQueue): Queue to hold list of digested reads.
         statq (multiprocessing.Queue, optional): Queue to use for aggregating statistics from digestion processes. Defaults to None.
         **digestion_kwargs: Kwargs passed to DigestedRead.
        
        Raises:
            KeyError: digestion_kwargs must contain: cutsite
        """    

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
        """
        Performs read digestion.
        
        Reads to digest are pulled from inq, digested with the DigestedRead class
        and the results placed on outq for writing.

        If a statq is provided, read digestion stats are placed into this queue for 
        aggregation.

        """        
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

  