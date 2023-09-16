#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith

Script generates a bed file of restriction fragment locations in a given genome.

"""
import pysam
import xopen
from typing import Iterator
import os
from loguru import logger
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
    **kwargs,
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

    from capcruncher_tools.api import digest_genome
    import polars as pl

    logger.info("Digesting genome")
    digest_genome(
        fasta=input_fasta,
        output=output_file,
        restriction_enzyme=recognition_site,
        remove_recognition_site=True,
        minimum_slice_length=18,
        n_threads=1,
    )

    logger.info("Digestion complete")

    if sort:
        logger.info("Sorting output")
        df = pl.read_csv(
            output_file, separator="\t", new_columns=["chrom", "start", "end", "name"]
        )

        # If changing the order, also need to change the fragment number
        df = (
            df.sort(["chrom", "start"])
            .drop(columns="name")
            .with_row_count("name")[["chrom", "start", "end", "name"]]
        )

        df.write_csv(output_file, separator="\t", has_header=False)
