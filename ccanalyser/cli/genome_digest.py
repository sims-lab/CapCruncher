#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith

Script generates a bed file of restriction fragment locations in a given genome.

"""
import click
import pysam
import xopen
from ccanalyser.cli import cli
from ccanalyser.tools.digest import DigestedChrom
from ccanalyser.utils import get_re_site


def parse_chromosomes(fasta: pysam.FastxFile) -> pysam.FastqProxy:

    for chrom in pysam.FastxFile(fasta):
        yield chrom


@cli.command()
@click.argument('input_fasta')
@click.option('-r', '--recognition_site', required=True)
@click.option('-l', '--logfile', default='genome_digest.log')
@click.option('-o', '--output_file', default='genome_digested.bed')
@click.option('--remove_cutsite', default=True)
def genome_digest(
    input_fasta,
    recognition_site=None,
    logfile=None,
    output_file=None,
    remove_cutsite=True,
):

     
    '''Digests the supplied genome fasta file and generates a bed file containing the locations of all restriction fragments
       produced by the supplied restriction enzyme '''
    #TODO: Include option to keep or remove the cutsite. For now will just remove to keep inline with the fastq digestion script

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
    
    with xopen(logfile, 'w') as output:
        for chrom, n_fragments in fragment_stats.items():
            output.write(f'{chrom}\t{n_fragments}\n')
