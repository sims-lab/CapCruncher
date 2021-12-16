#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith

Script generates a bed file of restriction fragment locations in a given genome.

"""
import pysam
import xopen
from capcruncher.tools.digest import DigestedChrom
from capcruncher.utils import get_re_site
from typing import Iterator
import os
import logging
import pandas as pd


def parse_chromosomes(fasta: pysam.FastxFile) -> Iterator[pysam.FastqProxy]:
    """Parses a whole genome fasta file and yields chromosome entries.

    Args:
        fasta (pysam.FastxFile): Fasta file to process.

    Yields:
        Iterator[pysam.FastqProxy]: Chromosome entry.
    """

    for chrom in pysam.FastxFile(fasta):
        yield chrom


def digest(
    input_fasta: os.PathLike,
    recognition_site: str,
    logfile: os.PathLike = "genome_digest.log",
    output_file: os.PathLike = "genome_digest.bed",
    remove_cutsite: bool = True,
    sort=False,
):
    """
    Performs in silico digestion of a genome in fasta format.

    Digests the supplied genome fasta file and generates a bed file containing the
    locations of all restriction fragments produced by the supplied restriction enzyme.

    A log file recording the number of restriction fragments for the suplied genome is also
    generated.

    \f
    Args:
     input_fasta (os.PathLike): Path to fasta file containing whole genome sequence, split by chromosome
     recognition_site (str): Restriction enzyme name/ Sequence of recognition site.
     logfile (os.PathLike, optional): Output path of the digestion logfile. Defaults to genome_digest.log.
     output_file (os.PathLike, optional): Output path for digested chromosome bed file. Defaults to genome_digest.bed.
     remove_cutsite (bool, optional): Determines if restriction site is removed. Defaults to True.
    """

    # TODO: Include option to keep or remove the cutsite. For now will just remove to keep inline with the fastq digestion script

    fragment_stats = dict()
    fragment_number = 0
    cut_sequence = get_re_site(recognition_site=recognition_site)

    with xopen.xopen(output_file, "w") as output:

        for chrom in parse_chromosomes(input_fasta):

            logging.info(f"Processing chrom {chrom.name}")

            digested_chrom = DigestedChrom(
                chrom,
                cut_sequence,
                fragment_number_offset=fragment_number,
                fragment_min_len=1,
            )

            for n_fragments, fragment in enumerate(digested_chrom.fragments):
                if n_fragments % 10000 == 0:
                    logging.info(f"Written {n_fragments} fragments")

                output.write(fragment)

            fragment_stats[chrom.name] = n_fragments
            fragment_number += n_fragments + 1

    if sort:
        logging.info("Sorting output")
        df = pd.read_csv(output_file, sep="\t", names=["chrom", "start", "end", "name"])

        # If changing the order, also need to change the fragment number
        df = (
            df.sort_values(["chrom", "start"])
            .drop(columns="name")
            .reset_index(drop=True)
            .reset_index()
            .rename(columns={"index": "name"})[["chrom", "start", "end", "name"]]
        )

        df.to_csv(output_file, sep="\t", index=False, header=False)

    with xopen.xopen(logfile, "w") as output:
        for chrom, n_fragments in fragment_stats.items():
            output.write(f"{chrom}\t{n_fragments}\n")
