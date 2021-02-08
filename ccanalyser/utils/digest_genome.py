#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith

Script generates a bed file of restriction fragment locations in a given genome.

"""
import re
import pysam
from ccanalyser.utils.helpers import get_re_site
import xopen


def parse_chromosomes(fasta: pysam.FastxFile) -> pysam.FastqProxy:

    for chrom in pysam.FastxFile(fasta):
        yield chrom


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


def main(
    input_fasta,
    recognition_site=None,
    logfile=None,
    output_file=None,
    remove_cutsite=True,
):

     
    
    #TODO: Include option to keep or remove the cutsite. For now will just remove to keep inline with the fastq digestion script


    # if not re.match(r'[GgAaTtCc]+', recognition_site):
    #     re_fragment_name = recognition_site
    # else:
    #     re_fragment_name = "U"

    # with open(logfile, "w") as log, open(output_file, "w") as bed_out, FastxFile(
    #     input_fasta
    # ) as fasta_file:


    fragment_stats = dict()
    fragment_number_offset = 0
    cut_sequence = get_re_site(recognition_site=recognition_site)

    with xopen.xopen(output_file, 'w') as output:
        
        for chrom in parse_chromosomes(input_fasta):

            digested_chrom = DigestedChrom(chrom, 
                                           cut_sequence, 
                                           fragment_number_offset=fragment_number_offset,
                                           fragment_min_len=1)
        
            for n_fragments, fragment in enumerate(digested_chrom.fragments):

                output.write(fragment)
            
            fragment_stats[chrom.name] = n_fragments
            fragment_number_offset += n_fragments + 1
    
    with open(logfile, 'w') as output:
        for chrom, n_fragments in fragment_stats.items():
            output.write(f'{chrom}\t{n_fragments}\n')

