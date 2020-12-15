#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith

Script generates a bed file of restriction fragment locations in a given genome.

"""
import re
from pysam import FastxFile
from ccanalyser.utils.helpers import get_re_site


def main(
    input_fasta,
    recognition_site=None,
    logfile=None,
    output_file=None,
    remove_cutsite=True,
):

    
    cut_sequence = get_re_site(recognition_site=recognition_site)
    cut_sequence_len = len(cut_sequence)
    #TODO: Include option to keep or remove the cutsite. For now will just remove to keep inline with the fastq digestion script



    if not re.match(r'[GgAaTtCc]+', recognition_site):
        re_fragment_name = recognition_site
    else:
        re_fragment_name = "U"

    with open(logfile, "w") as log, open(output_file, "w") as bed_out, FastxFile(
        input_fasta
    ) as fasta_file:

        for seq_entry in fasta_file:

            # Find match positions of the restriction enzyme sequence
            seq_length = len(seq_entry.sequence)
            match_positions = [
                m.start() for m in re.finditer(cut_sequence, seq_entry.sequence.upper())
            ]

            # iterate through matches and write to bed file
            slice_start = 0
            for match_index, match_pos in enumerate(match_positions):
                if remove_cutsite and not slice_start == 0:
                    slice_start += cut_sequence_len
                
                slice_end = match_pos

                if slice_start != slice_end:
                    slice_name = f"{re_fragment_name}_{seq_entry.name}_{match_index}"
                    bed_out.write(
                        f"{seq_entry.name}\t{slice_start}\t{slice_end}\t{slice_name}\n"
                    )
                slice_start = slice_end

            # handle last slice
            if slice_start != seq_length:
                slice_end = seq_length
                slice_name = f"{re_fragment_name}_{seq_entry.name}_{match_index}"
                bed_out.write(
                    f"{seq_entry.name}\t{slice_start}\t{slice_end}\t{slice_name}\n"
                )

            # Print total slices per chr to log file
            log.write(f"{seq_entry.name}: {len(match_positions)}\n")
