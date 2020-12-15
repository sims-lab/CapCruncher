import multiprocessing as mp
import re
from collections import Counter

import pandas as pd
import numpy as np
from pysam import FastxFile
from xopen import xopen
from ccanalyser.utils.helpers import get_re_site
from random import randint
from multiprocessing import Queue, Pool, Manager, Process
from typing import Union



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
        self.slices_total_counter = len(self.slice_indexes) - 1
        self.slices_valid_counter = 0
        self.has_slices = self.slices_total_counter > 1
        self.slices = self.get_slices()

    def get_recognition_site_indexes(self):
        indexes = [
            re_site.start()
            for re_site in self.recognition_re.finditer(self.read.sequence.upper())
        ]

        indexes.insert(0, 0)
        indexes.append(len(self.read.sequence))

        return indexes

    def get_slices(self):

        indexes = self.slice_indexes
        slice_no = self.slice_number_offset
        slices_list = []

        if self.has_slices or self.allow_undigested:

            # Iterate through offset indexes to get correct start and end
            for ii, (slice_start, slice_end) in enumerate(zip(indexes, indexes[1:])):

                # If this is not the first slice
                if ii > 0:
                    slice_start += self.recognition_len

                if self.is_valid_slice(slice_start, slice_end):
                    slices_list.append(
                        self.prepare_slice(slice_start, slice_end, slice_no)
                    )

                    self.slices_valid_counter += 1
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

    def is_valid_slice(self, start, end):
        if (end - start) >= self.min_slice_length:
            return True

    def __str__(self):
        return ("\n".join(self.slices) + "\n") if self.slices else ""



