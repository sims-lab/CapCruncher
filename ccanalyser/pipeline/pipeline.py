#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Capture-C Pipeline
==================

This script processes data from the capture-c or NG capture-c sequencing
protocols designed to identify 3d interactions in the genome.

It takes Illumina paired-end sequencing reads in fastq format (gz compressed)
as input and performs the following steps:

1: Identifies all restriction fragments in the genome (digest_genome.py)
2: Quality control of raw reads (fastqc, multiqc)
3: Removal of PCR duplicates (based on exact sequence matches; deduplicate_fastq.py)
4: Trimming of reads to remove adaptor sequence (trim_galore)
5: Combining overlapping read pairs (flash)
6: Splits fastq(s) into smaller files (enables fast parallel processing; split_fastq.py)
7: In silico digestion of reads (digest_fastq.py)
8: Alignment of fastq files with a user specified aligner (i.e. bowtie/bowtie2; BWA is not supported)
9: Analysis of alignment statistics (picard CollectAlignmentSummaryMetrics, multiqc)
10: Intersection of mapped reads with: capture probes, exclusion regions, blacklist, restriction fragments
11: Classification of mapped read slices (ccanalyser.py)
12: Generation of bedgraphs/BigWigs (convert_tsv_to_bedgraph.py and make_bigwig)
13: Generation of a UCSC track hub for visualisation
14: Collation of run statistics and generation of a run report


@authors: asmith, dsims
"""

from collections import defaultdict
import itertools
import os
import re
import sys
import glob
from cgatcore.pipeline.parameters import PARAMS
from ruffus.task import originate

# import seaborn as sns
import trackhub
from cgatcore import pipeline as P
from cgatcore.iotools import zap_file
from pybedtools.helpers import get_chromsizes_from_ucsc
import pandas as pd
from ruffus import (
    active_if,
    add_inputs,
    collate,
    follows,
    merge,
    mkdir,
    regex,
    transform,
    suffix,
)
from ccanalyser.tools.statistics import (
    collate_slice_data,
    collate_read_data,
    collate_cis_trans_data,
    collate_histogram_data,
    extract_trimming_stats,
)
import click
from ccanalyser.cli import cli
from ccanalyser.tools.pileup import CCBedgraph
from ccanalyser.utils import is_on, is_none, zap_files

##############################
#   Set-up global parameters #
##############################

P.get_parameters("config.yml")

##############################
#  Pipeline stages           #
##############################


@follows(mkdir("pre_ccanalysis"), mkdir("pre_ccanalysis/fastqc"))
@transform("*.fastq*", regex(r"(.*).fastq.*"), r"pre_ccanalysis/fastqc/\1_fastqc.zip")
def qc_reads(infile, outfile):
    """Quality control of raw sequencing reads"""
    outdir = os.path.dirname(outfile)
    statement = " ".join(
        [
            "fastqc",
            "%(infile)s",
            "-q",
            "-t",
            "%(pipeline_n_cores)s",
            "--nogroup",
            "--outdir",
            "%(outdir)s",
        ]
    )

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("statistics"))
@merge(qc_reads, "statistics/fastqc_report.html")
def multiqc_reads(infile, outfile):
    """Collate fastqc reports into single report using multiqc"""

    bn = os.path.basename(outfile)
    dn = os.path.dirname(outfile)

    s1 = "rm -f %(outfile)s &&"
    s2 = "export LC_ALL=en_US.UTF-8 &&"
    s3 = "export LANG=en_US.UTF-8 &&"
    s4 = " ".join(["multiqc", "pre_ccanalysis/fastqc/", "-o", "%(dn)s", "-n", "%(bn)s"])

    statement = " ".join([s1, s2, s3, s4])

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("pre_ccanalysis/restriction_enzyme_map/"))
@transform(
    P.PARAMS.get("genome_fasta"),
    regex(r".*/(.*).fa.*"),
    r"pre_ccanalysis/restriction_enzyme_map/genome.digest.bed.gz",
)
def digest_genome(infile, outfile):
    """Digest genome using restriction enzyme and output fragments in bed file"""

    assert infile, "genome_fasta not provided, please provide path in config.yml"

    tmp = outfile.replace('.gz', '')
    statement = """ccanalyser genome-digest
                   %(infile)s
                   -l %(tmp)s.log
                   -o %(tmp)s
                   -r %(analysis_restriction_enzyme)s
                   && cat %(tmp)s | sort -k1,1 -k2,2n > %(tmp)s.sorted
                   && mv %(tmp)s.sorted %(tmp)s
                   && pigz -p 4 %(tmp)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("pre_ccanalysis/split"))
@collate(
    "*.fastq.gz",
    regex(r"(.*)_R*[12].fastq.*"),
    r"pre_ccanalysis/split/\1.log",
)
def split_fastq(infiles, outfile):
    """Splits the combined (flashed) fastq files into chunks for parallel processing"""

    infiles = " ".join(infiles)
    output_prefix = outfile.replace(".log", "")

    statement = """ccanalyser fastq-split
                %(infiles)s
                -m unix
                -o %(output_prefix)s
                -n %(split_n_reads)s
                --no-gzip
                > %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(is_on(P.PARAMS.get("deduplication_pre-dedup")))
@follows(
    mkdir("pre_ccanalysis/deduplicated"),
    mkdir("pre_ccanalysis/deduplicated/deduplicated_ids"),
    split_fastq,
)
@collate(
    "pre_ccanalysis/split/*.fastq*",
    regex(r"pre_ccanalysis/split/(.*)_part(\d+)_[12].fastq(?:.gz)?"),
    r"pre_ccanalysis/deduplicated/deduplicated_ids/\1_\2.json.gz",
    extras=[r'\1', r'\2']
)
def parse_duplicate_reads(infiles, outfile, sample_name, part_no):

    """Checks for duplicate read1/read2 pairs in a pair of fastq files
    any duplicates are discarded"""

    fq1, fq2 = [os.path.abspath(fn) for fn in infiles]

    statement = """ccanalyser fastq-deduplicate parse
                   %(fq1)s %(fq2)s
                   -o %(outfile)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="6G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@collate(
    parse_duplicate_reads,
    regex(r"pre_ccanalysis/deduplicated/deduplicated_ids/(.*)_\d*.json.gz"),
    r"pre_ccanalysis/deduplicated/deduplicated_ids/\1.json.gz",
)
def identify_duplicate_reads(infiles, outfile):
    
    infiles_str = " ".join(infiles)
    statement = """ccanalyser fastq-deduplicate identify
                   %(infiles_str)s
                   -o %(outfile)s
                 """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="32G",
        job_condaenv=P.PARAMS["conda_env"],
    )

    for fn in infiles:
        zap_file(fn)

@active_if(is_on(P.PARAMS.get("deduplication_pre-dedup")))
@follows(
    parse_duplicate_reads,
    identify_duplicate_reads,
    mkdir("statistics/deduplication/data/"),
)
@collate(
    "pre_ccanalysis/split/*.fastq*",
    regex(r".*/(.*_part\d+)_[12].fastq(?:.gz)?"),
    r"pre_ccanalysis/deduplicated/\1_1.fastq",
)
def remove_duplicate_reads(infiles, outfile):

    """Checks for duplicate read1/read2 pairs in a pair of fastq files
    any duplicates are discarded"""

    fq1, fq2 = infiles
    sample = re.match(r".*/(.*)(_part\d+)_[12].fastq(?:.gz)?", fq1)
    sample_name = sample.group(1)
    sample_part = sample.group(2)
    dd_ids = f"pre_ccanalysis/deduplicated/deduplicated_ids/{sample_name}.json.gz"
    stats_prefix = f"statistics/deduplication/data/{sample_name}_{sample_part}"
    output_prefix = outfile.replace("_1.fastq", "")

    statement = """ccanalyser fastq-deduplicate remove
                            %(fq1)s %(fq2)s
                            -d %(dd_ids)s
                            -o %(output_prefix)s
                            --sample_name %(sample_name)s
                            --stats_prefix %(stats_prefix)s
                            """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="6G",
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Replace infiles with empty files
    for fn in infiles:
        zap_file(fn)



@follows(remove_duplicate_reads)
@merge(
    "statistics/deduplication/data/*.csv",
    "statistics/deduplication/deduplication.reads.csv",
)
def collate_deduplication_stats(infiles, outfile):

    stats_prefix = outfile.replace(".reads.csv", "")

    df_stats = collate_read_data(infiles)

    df_stats_read = df_stats.query('stat_type != "reads_removed"')

    df_stats.to_csv(f"{stats_prefix}.summary.csv", index=False)
    df_stats_read.to_csv(
        outfile, index=False
    )  # Modified to enable more streamlined summary at final stage


@follows(
    mkdir("pre_ccanalysis/trimmed"),
    remove_duplicate_reads,
    mkdir("statistics/trimming/data/"),
)
@collate(
    "pre_ccanalysis/deduplicated/*.fastq*",
    regex(r"pre_ccanalysis/deduplicated/(.*)_[12].fastq(?:.gz)?"),
    r"pre_ccanalysis/trimmed/\1_1_val_1.fq.gz",
)
def trim_reads(infiles, outfile):
    """Trim adaptor sequences using Trim-galore"""

    fq1, fq2 = infiles
    fq1_basename, fq2_basename = os.path.basename(fq1), os.path.basename(fq2)

    outdir = os.path.dirname(outfile)
    trim_options = (
        P.PARAMS["trim_options"] if not is_none(P.PARAMS["trim_options"]) else ""
    )
    statement = """trim_galore
                   --cores %(pipeline_n_cores)s
                   --paired %(trim_options)s
                   --gzip
                   -o %(outdir)s
                   %(fq1)s
                   %(fq2)s
                   && mv pre_ccanalysis/trimmed/%(fq1_basename)s_trimming_report.txt statistics/trimming/data
                   && mv pre_ccanalysis/trimmed/%(fq2_basename)s_trimming_report.txt statistics/trimming/data
                   """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


    for fn in infiles:
        zap_file(fn)

@follows(trim_reads)
@merge("statistics/trimming/data/*.txt", r"statistics/trimming/trimming.summary.csv")
def collate_trimming_statistics(infiles, outfile):

    trimming_stats = []
    for fn in infiles:
        stats = extract_trimming_stats(fn)
        trimming_stats.append(stats)

    df_trimming_stats = (
        pd.DataFrame(trimming_stats)
        .groupby(["sample", "read_number", "read_type"])
        .sum()
        .reset_index()
        .melt(
            id_vars=["sample", "read_number", "read_type"],
            var_name="stat_type",
            value_name="stat",
        )
    )

    df_trimming_stats.to_csv("statistics/trimming/trimming.summary.csv", index=False)


@follows(trim_reads, mkdir("pre_ccanalysis/flashed"))
@collate(
    "pre_ccanalysis/trimmed/*.fq*",
    regex(r"pre_ccanalysis/trimmed/(.*)_[12]_.*.fq(?:.gz)?"),
    r"pre_ccanalysis/flashed/\1.extendedFrags.fastq.gz",
)
def combine_reads(infiles, outfile):
    """Combine overlapping paired-end reads using flash"""
    fq1, fq2 = infiles
    output_prefix = outfile.replace(".extendedFrags.fastq.gz", "")
    statement = """flash
                   -z
                   -t %(pipeline_n_cores)s
                   -o %(output_prefix)s
                    %(fq1)s
                    %(fq2)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(
    mkdir("pre_ccanalysis/digested"), combine_reads, mkdir("statistics/digestion/data")
)
@transform(
    "pre_ccanalysis/flashed/*.fastq.gz",
    regex(r"pre_ccanalysis/flashed/(.*).extendedFrags.fastq.gz"),
    r"pre_ccanalysis/digested/\1.flashed.fastq.gz",
)
def digest_flashed_reads(infile, outfile):
    """In silico restriction enzyme digest of combined (flashed) read pairs"""

    fn = os.path.basename(infile).replace(".extendedFrags.fastq.gz", "")
    statement = """ccanalyser fastq-digest
                   %(infile)s
                   -m flashed
                   -r %(analysis_restriction_enzyme)s
                   -o %(outfile)s
                   --stats_prefix %(outfile)s.stats
                   --minimum_slice_length 18
                   --compression_level %(pipeline_advanced_compression)s
                   -p 1
                   --stats_prefix statistics/digestion/data/%(fn)s.flashed
                    """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=4,
        job_condaenv=P.PARAMS["conda_env"],
    )

    zap_file(infile)



@follows(combine_reads, mkdir("statistics/digestion"))
@collate(
    "pre_ccanalysis/flashed/*.fastq.gz",
    regex(r"pre_ccanalysis/flashed/(.*).notCombined_[12].fastq.gz"),
    r"pre_ccanalysis/digested/\1.pe.fastq.gz",
)
def digest_pe_reads(infiles, outfile):
    """In silico restriction enzyme digest of non-combined (non-flashed) read pairs"""

    fq1, fq2 = infiles
    fn = os.path.basename(fq1).replace(".notCombined_1.fastq.gz", "")

    statement = """ccanalyser fastq-digest
                   %(fq1)s
                   %(fq2)s
                   -m pe
                   -r %(analysis_restriction_enzyme)s
                   -o %(outfile)s
                   --minimum_slice_length 18
                   -p 1
                   --compression_level %(pipeline_advanced_compression)s
                   --stats_prefix statistics/digestion/data/%(fn)s.pe
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=4,
        job_condaenv=P.PARAMS["conda_env"],
    )

    for fn in infiles:
        zap_file(fn)

@follows(digest_flashed_reads, digest_pe_reads)
@merge("statistics/digestion/data/*", "statistics/digestion/digestion.reads.csv")
def collate_digestion_stats(infiles, outfile):

    stats_prefix = outfile.replace(".reads.csv", "")
    data = defaultdict(list)

    for fn in infiles:
        if ".filtered.histogram" in fn:
            data["hist_filt"].append(fn)
        elif ".unfiltered.histogram" in fn:
            data["hist_unfilt"].append(fn)
        elif ".slice" in fn:
            data["slice"].append(fn)
        elif ".read" in fn:
            data["read"].append(fn)

    df_hist_filt = collate_histogram_data(data["hist_filt"])
    df_hist_unfilt = collate_histogram_data(data["hist_unfilt"])
    df_slice = collate_read_data(data["slice"])
    df_read = collate_read_data(data["read"])

    df_hist = pd.concat(
        [df_hist_unfilt.assign(filtered=0), df_hist_filt.assign(filtered=1)]
    ).sort_values(["sample", "read_type", "number_of_slices"])

    df_hist.to_csv(f"{stats_prefix}.histogram.csv", index=False)
    df_slice.to_csv(f"{stats_prefix}.slice.csv", index=False)
    df_read.to_csv(outfile, index=False)


@follows(digest_flashed_reads, digest_pe_reads)
def fastq_preprocessing():
    pass


@follows(mkdir("pre_ccanalysis"), fastq_preprocessing)
@transform(
    [digest_flashed_reads, digest_pe_reads],
    regex(r"pre_ccanalysis/digested/(.*).fastq.gz"),
    r"pre_ccanalysis/aligned/\1.bam",
)
def align_slices(infile, outfile):
    """ Aligns digested fq files using designated aligner. (Default: Bowtie2)"""

    aligner = P.PARAMS["align_aligner"]
    index_flag = (
        P.PARAMS["align_index_flag"]
        if not is_none(P.PARAMS["align_index_flag"])
        else ""
    )
    options = (
        P.PARAMS["align_options"] if not is_none(P.PARAMS["align_options"]) else ""
    )

    s1 = "%(aligner)s %(options)s %(index_flag)s %(genome_aligner_index)s %(infile)s"  # Align reads
    s2 = "| samtools view -b -S > %(outfile)s"  # Convert to bam
    s3 = "&& samtools sort %(outfile)s -o %(outfile)s.sorted.bam -m 2G -@ %(pipeline_n_cores)s"  # Sort reads
    s4 = "&& mv -f %(outfile)s.sorted.bam %(outfile)s"  # Replace unsorted reads with sorted

    statement = " ".join([s1, s2, s3, s4])

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory="4G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@collate(
    align_slices,
    regex(r"pre_ccanalysis/aligned/(.*)_part\d+.*.bam"),
    r"pre_ccanalysis/aligned/\1.bam",
)
def merge_bam_files(infiles, outfile):
    """Combines bam files (by flashed/non-flashed status and sample)"""
    fnames = " ".join(infiles)

    statement = """samtools merge %(outfile)s %(fnames)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform([align_slices, merge_bam_files], regex("(.*).bam"), r"\1.bam.bai")
def index_bam(infile, outfile):
    statement = """samtools index %(infile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="1G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("statistics/mapping_statistics"), index_bam, align_slices)
@transform(
    merge_bam_files,
    regex(r"pre_ccanalysis/(.*).bam"),
    r"statistics/mapping_statistics/\1.picard.metrics",
)
def mapping_qc(infile, outfile):
    """Uses picard CollectAlignmentSummaryMetrics to get mapping information."""

    cmd = [
        "picard",
        "CollectAlignmentSummaryMetrics",
        "VALIDATION_STRINGENCY=LENIENT",
        "R=%(genome_fasta)s",
        "I=%(infile)s",
        "O=%(outfile)s",
        "&> %(outfile)s.log",
    ]

    statement = " ".join(cmd)

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@merge(mapping_qc, "statistics/mapping_report.html")
def mapping_multiqc(infiles, outfile):
    """Combines mapping metrics using multiqc"""

    indir = os.path.dirname(infiles[0])
    out_fn = os.path.basename(outfile)
    out_dn = os.path.dirname(outfile)
    statement = """rm -f %(outfile)s &&
                   export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc
                   %(indir)s
                   -o %(out_dn)s
                   -n %(out_fn)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="8G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("ccanalysis/annotations"))
@transform(
    align_slices,
    regex(r"pre_ccanalysis/aligned/(.*).bam"),
    r"ccanalysis/annotations/\1.bam.bed",
)
def bam_to_bed(infile, outfile):
    """Converts bam files to bed for faster intersection"""


    statement = """bedtools bamtobed -i %(infile)s | sort -k1,1 -k2,2n > %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(
    P.PARAMS.get("analysis_capture_oligos"),
    regex(r".*"),
    r"ccanalysis/annotations/exclude.bed",
)
def build_exclusion_bed(infile, outfile):
    """Generates exclusion window around each capture site"""

    statement = """bedtools slop
                    -i %(infile)s -g %(genome_chrom_sizes)s -b %(analysis_reporter_exclusion_zone)s
                    | bedtools subtract -a - -b %(analysis_capture_oligos)s
                    | sort -k1,1 -k2,2n 
                    > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@originate("ccanalysis/annotations/capture.bed")
def sort_capture_oligos( outfile):
    
    """Sorts the capture oligos for bedtools intersect with --sorted option"""

    statement = """cat %(analysis_capture_oligos)s | sort -k1,1 -k2,2n > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@active_if(P.PARAMS.get('analysis_blacklist'))
@originate("ccanalysis/annotations/capture.bed")
def sort_blacklist(outfile):
    
    """Sorts the capture oligos for bedtools intersect with --sorted option"""

    statement = """cat %(analysis_blacklist)s | sort -k1,1 -k2,2n > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(digest_genome, build_exclusion_bed, sort_capture_oligos, sort_blacklist)
@transform(
    bam_to_bed,
    regex(r"ccanalysis/annotations/(.*).bam.bed"),
    add_inputs(
        [
            {
                "name": "restriction_fragment",
                "fn": "pre_ccanalysis/restriction_enzyme_map/genome.digest.bed.gz",
                "action": "get",
                "fraction": 0.2,
            },
            {
                "name": "capture",
                "fn": 'ccanalysis/annotations/capture.bed',
                "action": "get",
                "fraction": 0.9,
            },
            {
                "name": "exclusion",
                "fn": "ccanalysis/annotations/exclude.bed",
                "action": "get",
                "fraction": 1e-9,
            },
            {
                "name": "exclusion_count",
                "fn": "ccanalysis/annotations/exclude.bed",
                "action": "count",
                "fraction": 1e-9,
            },
            {
                "name": "capture_count",
                "fn": 'ccanalysis/annotations/capture.bed',
                "action": "count",
                "fraction": 0.9,
            },
            {
                "name": "blacklist",
                "fn": 'ccanalysis/annotations/blacklist.bed',
                "action": "count",
                "fraction": 1e-9,
            },
        ]
    ),
    r"ccanalysis/annotations/\1.annotations.tsv",
)
def annotate_slices(infile, outfile):

    slices = infile[0]
    flags = {"name": "-n", "fn": "-b", "action": "-a", "fraction": "-f"}

    cmd_args = []
    for args in infile[1]:
        for arg_name, arg in args.items():
            cmd_args.append(f'{flags.get(arg_name)} {arg if arg else "-"}')

    cmd_args = " ".join(cmd_args)
    statement = """ccanalyser slices-annotate
                    %(slices)s
                    %(cmd_args)s
                    -o %(outfile)s
                """

    P.run(
        statement.replace("\n", " "),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=4,
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zero bed file used
    zap_file(infile[0])


@follows(fastq_preprocessing, annotate_slices)
def pre_ccanalysis():
    pass


@follows(
    pre_ccanalysis,
    annotate_slices,
    mkdir("ccanalysis/reporters/unfiltered/"),
    mkdir("statistics/ccanalysis/data"),
)
@transform(
    align_slices,
    regex(r"pre_ccanalysis/aligned/(.*).bam"),
    add_inputs(r"ccanalysis/annotations/\1.annotations.tsv"),
    r"ccanalysis/reporters/unfiltered/\1.log",
)
def ccanalyser(infiles, outfile):
    """Processes bam files and annotations, filteres slices and outputs
    reporter slices for each capture site"""

    bam, annotations = infiles
    sample = re.match(r".*/(.*)_(part\d+).(flashed|pe).bam", bam)
    sample_name = sample.group(1)
    sample_part = sample.group(2)
    sample_read_type = sample.group(3)

    output_prefix = outfile.replace(".log", "")
    stats_prefix = (
        f"statistics/ccanalysis/data/{sample_name}_{sample_part}_{sample_read_type}"
    )

    statement = """ccanalyser reporters-identify
                   %(analysis_method)s
                   -b %(bam)s
                   -a %(annotations)s
                   -o %(output_prefix)s
                   --stats_prefix %(stats_prefix)s
                   --sample_name %(sample_name)s
                   --read_type %(sample_read_type)s
                   > %(outfile)s 2>&1"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_advanced_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zero annotations
    zap_file(annotations)



@follows(mkdir("ccanalysis/reporters/combined"), ccanalyser)
@collate(
    "ccanalysis/reporters/unfiltered/*.tsv",
    regex(
        r".*/(?P<sample>.*)_part\d+.(flashed|pe).(?P<capture>.*).(slices|fragments).tsv"
    ),
    r"ccanalysis/reporters/combined/\1.\2.\3.\4.tsv",
    extras=[r"\1", r"\2", r"\3", r'\4'],
)
def collate_reporters(infiles, outfile, *grouping_args):

    # Need to concat tsv files but remove headers, the sed command performs the removal
    statement = f"cat {' '.join(infiles)} | sed -e '1n' -e '/.*parent_read.*/d' > {outfile}"

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zero un-aggregated reporters
    for fn in infiles:
        zap_file(fn)

@follows(ccanalyser, mkdir("ccanalysis/reporters/deduplicated"))
@transform(collate_reporters, 
           regex(r".*/(?P<sample>.*).(flashed|pe).(?P<capture>.*).fragments.tsv"),
           r'ccanalysis/reporters/deduplicated/\1.\2.\3.json.gz',
           extras=[r'\2']
)
def deduplicate_fragments(infile, outfile, read_type):

    statement = """ccanalyser reporters-deduplicate
                   identify
                   %(infile)s
                   --read_type %(read_type)s
                   -o %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory="32G",
        job_condaenv=P.PARAMS["conda_env"],
    )

@follows(deduplicate_fragments)
@transform(
    collate_reporters,
    regex(r"ccanalysis/reporters/combined/(.*)\.(flashed|pe)\.(.*)\.slices.tsv"),
    add_inputs(r"ccanalysis/reporters/deduplicated/\1.\2.\3.json.gz"),
    r"ccanalysis/reporters/deduplicated/\1.\2.\3.slices.tsv",
    extras=[r'\1', r'\2', r'\3']
)
def deduplicate_slices(infile, outfile, sample_name, read_type, capture_oligo):


    slices, duplicated_ids = infile
    stats_prefix = f"statistics/ccanalysis/data/{sample_name}_{read_type}_{capture_oligo}"

    statement = """ccanalyser reporters-deduplicate
                   remove
                   %(slices)s
                   -d %(duplicated_ids)s
                   -o %(outfile)s
                   --stats_prefix %(stats_prefix)s
                   --sample_name %(sample_name)s
                   --read_type %(read_type)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory=P.PARAMS["pipeline_advanced_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zero non-deduplicated reporters
    zap_file(slices)


@collate(
    deduplicate_slices,
    regex(
        r".*/(?P<sample>.*).(?:flashed|pe).(?P<capture>.*).slices.tsv"
    ),
    r"ccanalysis/reporters/\1.\2.tsv.gz",
    extras=[r"\1", r"\2"],
)
def collate_reporters_final(infiles, outfile, *grouping_args):

    # Need to concat tsv files but remove headers, the sed command performs the removal
    statement = f"cat {' '.join(infiles)} | sed -e '1n' -e '/.*parent_read.*/d' | pigz -p 4 > {outfile}"

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(collate_reporters_final)
@merge("statistics/ccanalysis/data/*", "statistics/ccanalysis/ccanalysis.reads.csv")
def collate_ccanalyser_stats(infiles, outfile):

    stats_prefix = outfile.replace(".reads.csv", "")
    data = defaultdict(list)

    for fn in infiles:
        if ".read.stats" in fn:
            data["read"].append(fn)
        elif ".slice.stats" in fn:
            data["slice"].append(fn)
        elif ".reporter.stats" in fn:
            data["reporter"].append(fn)

    # Slice data
    collate_slice_data(data["slice"]).to_csv(f"{stats_prefix}.slices.csv", index=False)

    # Reporter data
    collate_cis_trans_data(data["reporter"]).to_csv(
        f"{stats_prefix}.reporters.csv", index=False
    )

    # Read data
    collate_read_data(data["read"]).to_csv(outfile, index=False)


@follows(deduplicate_slices, collate_ccanalyser_stats)
def post_ccanalysis():
    pass


@follows(mkdir("ccanalysis/interactions"))
@transform(
    collate_reporters_final,
    regex(r"ccanalysis/reporters/(.*)\.(.*).tsv.gz"),
    r"ccanalysis/interactions/\1.\2.tsv.gz",
)
def count_interactions(infile, outfile):

    statement = [
        "ccanalyser interactions-count",
        '%(infile)s',
        "-o %(outfile)s",
        "--remove_exclusions",
        "--remove_capture" if P.PARAMS["analysis_method"] == "tri" else "",
        "> %(outfile)s.log",
    ]

    statement = " ".join(statement)

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=2,
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(
    count_interactions,
    regex(r"ccanalysis/interactions/(.*)\.(.*)\.tsv.gz"),
    add_inputs(digest_genome),
    r"ccanalysis/interactions/\1.\2.fragments.hdf5",
    extras=[r'\1', r'\2']
)
def store_interactions_at_fragment_level(infile, outfile, sample_name, capture_name):
    
    counts, rf_map = infile
    output_prefix = outfile.replace(f'.{capture_name}.fragments', '')

    statement = '''ccanalyser interactions-store
                   fragments
                   %(counts)s
                   -f %(rf_map)s
                   -g %(genome_name)s
                   -n %(capture_name)s
                   -o %(output_prefix)s
                   -c %(analysis_capture_oligos)s
                   --suffix fragments
                   '''
    
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@active_if(not is_none(P.PARAMS.get('plot_bin_size')))
@transform(
    store_interactions_at_fragment_level,
    regex(r"ccanalysis/interactions/(.*)\.(.*)\.fragments\.hdf5"),
    r"ccanalysis/interactions/\1.\2.log",
    extras=[r'\2']
)
def store_interactions_binned(infile, outfile, capture_name):


    bin_options = ' -b ' + " -b ".join(re.split(r"[,;]\s*|\s+", str(P.PARAMS["plot_bin_size"])))
    output_prefix = outfile.replace(f'.{capture_name}.log', '')

    statement =  f"""ccanalyser interactions-store
                  bins
                  %(infile)s
                  -o %(output_prefix)s
                  %(bin_options)s
                  --normalise
                  -p 8
                  > %(outfile)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@follows(store_interactions_at_fragment_level, store_interactions_binned)
@collate(
    'ccanalysis/interactions/*.hdf5',
    regex(r"ccanalysis/interactions/(.*)\.(.*)\.(?:fragments|\d+)\.hdf5"),
    r"ccanalysis/interactions/\1.hdf5",
    extras=[r'\1']
)
def merge_interactions(infiles, outfile, sample_name):

    infiles_str = ' '.join(infiles)
    statement =  f"""ccanalyser interactions-store
                     merge
                     %(infiles_str)s
                     -o %(outfile)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("ccanalysis/bedgraphs"))
@transform(merge_interactions, regex(r'.*/(.*).hdf5'), r'ccanalysis/bedgraphs/\1.raw.log', extras=[r'\1'])
def make_bedgraph_raw(infile, outfiles, sample_name):
    """Extract reporters in bedgraph format from stored interactions"""

    output_prefix = f'ccanalysis/bedgraphs/{sample_name}.raw'

    statement = '''ccanalyser interactions-bedgraph
                   %(infile)s
                   -o %(output_prefix)s
                   > %(output_prefix)s.log 
                   && pigz -p 8 %(output_prefix)s*.bedgraph'''
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@transform(merge_interactions, regex(r'.*/(.*).hdf5'), r'ccanalysis/bedgraphs/\1.normalised.log', extras=[r'\1'])
def make_bedgraph_normalised(infile, outfiles, sample_name):
    """Extract reporters in bedgraph format from stored interactions"""

    output_prefix = f'ccanalysis/bedgraphs/{sample_name}.normalised'

    statement = '''ccanalyser interactions-bedgraph
                   %(infile)s
                   -o %(output_prefix)s
                   --normalise
                   > %(output_prefix)s.log
                   && pigz -p 8 %(output_prefix)s*.bedgraph
                   '''
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@transform(merge_interactions, regex(r'.*/(.*).hdf5'), r'ccanalysis/bedgraphs/\1.windowed.log', extras=[r'\1'])
def make_bedgraph_windowed(infile, outfiles, sample_name):
    """Extract reporters in bedgraph format from stored interactions"""

    output_prefix = f'ccanalysis/bedgraphs/{sample_name}.windowed'

    statement = '''ccanalyser interactions-bedgraph
                   %(infile)s
                   -o %(output_prefix)s
                   --binsize 5000
                   --normalise
                   > %(output_prefix)s.log
                   && pigz -p 8 %(output_prefix)s*.bedgraph
                   '''
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )



@follows(mkdir("ccanalysis/bigwigs"), make_bedgraph_raw, make_bedgraph_normalised, make_bedgraph_windowed)
@transform(
    'ccanalysis/bedgraphs/*.gz',
    regex(r"ccanalysis/bedgraphs/(.*).bedgraph.gz"),
    r"ccanalysis/bigwigs/\1.bigWig",
)
def make_bigwig(infile, outfile):
    """Uses UCSC tools bedGraphToBigWig to generate bigWigs for each bedgraph"""

    tmp = infile.replace(".gz", "")


    #TODO: Bedgraphs are not being sorted according to unix conventions
    # Need to fix this but adding hack to sort this out
    statement = """zcat %(infile)s
                   | sort -k1,1 -k2,2n  > %(tmp)s
                   && bedGraphToBigWig %(tmp)s %(genome_chrom_sizes)s %(outfile)s
                   && rm %(tmp)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir('capture_compare/bedgraphs_union'))
@collate('ccanalysis/bedgraphs/*.gz', 
         regex(r'.*/(?:.*)\.(raw|normalised|windowed)\.(.*).bedgraph.gz'),
         r'capture_compare/bedgraphs_union/\2.\1.tsv',
         extras=[r'\1', r'\2'])
def make_union_bedgraph(infiles, outfile, normalisation_type, capture_name):

    infiles_str = ' '.join(infiles)
    sample_names = ' '.join([re.match(r'.*/(.*)\.(.*)\.(?:.*).bedgraph.gz', fn).group(1)
                            for fn in infiles])

    statement = '''bedtools 
                  unionbedg 
                  -i %(infiles_str)s 
                  -header 
                  -names %(sample_names)s
                  > %(outfile)s'''
    
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )














# @collate(
#     average_replicate_bedgraphs,
#     regex(r".*/(?:.*)\.(?P<capture>.*)\.mean\.bedgraph.gz"),
#     r"ccanalysis/bedgraph_compare/\1.log",
# )
# def compare_average_bedgraphs(infiles, outfile):

#     fn_regex = re.compile(r".*/(?P<sample>.*)\.(?P<capture>.*)\.mean\.bedgraph.gz")

#     if len(infiles) > 1:

#         for fn1, fn2 in itertools.combinations(infiles, 2):
#             bdg1, bdg2 = CCBedgraph(path=fn1), CCBedgraph(path=fn2)
#             bdg_sub = bdg1 - bdg2

#             sample1 = fn_regex.match(fn1)
#             sample2 = fn_regex.match(fn2)

#             out = f'{sample1.group("sample")}_vs_{sample2.group("sample")}.{sample1.group("capture")}.bedgraph.gz'
#             out_path = outfile.replace(sample1.group("capture") + ".log", out)
#             bdg_sub.to_file(out_path)











# @follows(store_interactions, mkdir("visualise"), mkdir("visualise/interaction_plots"))
# @transform(
#     "ccanalysis/interaction_counts_combined/*.hdf5",
#     regex(r"ccanalysis/interaction_counts_combined/(.*)_(.*)_binned.hdf5"),
#     add_inputs(r"ccanalysis/interaction_counts_combined/\1_restriction_fragments.hdf5"),
#     r"visualise/interaction_plots/\1.log",
# )
# def plot_interactions(infile, outfile):

#     if not is_none(P.PARAMS.get("plot_coordinates")):

#         cooler_binned, cooler_rf = infile
#         output_prefix = outfile.replace(".log", "")
#         normalisation = (
#             P.PARAMS.get("plot_advanced_normalisation")
#             if not is_none(P.PARAMS.get("plot_advanced_normalisation"))
#             else "infer"
#         )
#         normalisation_options = (
#             P.PARAMS.get("plot_advanced_normalisation_options")
#             if not is_none(P.PARAMS.get("plot_advanced_normalisation_options"))
#             else ""
#         )

#         statement = f"""ccanalyser ccanalysis plot_interactions
#                        {cooler_rf}
#                        {cooler_binned}
#                        -c %(plot_coordinates)s
#                        --method %(analysis_method)s
#                        --normalisation %(normalisation)s
#                        %(normalisation_options)s
#                        --cmap {P.PARAMS.get('plot_cmap', 'jet')}
#                        --thresh {P.PARAMS.get('plot_thresh', '0')}
#                        -o %(output_prefix)s
#                        -f %(plot_format)s
#                        > %(outfile)s 2>&1"""

#         P.run(
#             statement,
#             job_queue=P.PARAMS["pipeline_cluster_queue"],
#             job_condaenv=P.PARAMS["conda_env"],
#         )

#     else:
#         print("Not plotting as no coordinates provided")


@merge(
    [collate_deduplication_stats, collate_digestion_stats, collate_ccanalyser_stats],
    "statistics/run_statistics.csv",
)
def merge_run_statistics(infiles, outfile):

    df = pd.concat([pd.read_csv(fn) for fn in infiles])
    df.sort_values(
        ["sample", "read_type", "stat"], ascending=[True, True, False]
    ).to_csv(outfile)


# @follows(collate_deduplication_stats, collate_digestion_stats, collate_ccanalyser_stats)
# @merge(
#     [
#         "statistics/deduplication/*.csv",
#         "statistics/digestion/*.csv",
#         "statistics/ccanalysis/*.csv",
#     ],
#     "statistics/visualise_statistics.html",
# )
# def build_report(infile, outfile):
#     """Run jupyter notebook for reporting and plotting. First moves the notebook
#     then converts to html"""

#     statement = """rm statistics/visualise_statistics* -f &&
#                    papermill
#                    %(PACKAGE_DIR)s/stats/visualise_capture-c_stats.ipynb
#                    statistics/visualise_statistics.ipynb
#                    -p directory $(pwd)/statistics/ &&
#                    jupyter nbconvert
#                    --no-input
#                    --to html
#                    statistics/visualise_statistics.ipynb
#                    statistics/visualise_statistics.html
#                    """

#     P.run(
#         statement,
#         job_queue=P.PARAMS["pipeline_cluster_queue"],
#         job_condaenv=P.PARAMS["conda_env"],
#     )


# @active_if(is_on(P.PARAMS["hub_create"]))
# @merge([make_bigwig, build_report], os.path.join(P.PARAMS["hub_dir"], "hub.txt"))
# def create_trackhub(infiles, outfile):
#     def get_ucsc_color(color):
#         return ",".join([f"{i * 255}" for i in color])

#     def get_colors(items, colors=None):

#         if not colors:
#             colors = sns.color_palette("rainbow", len(items))
#             return [get_ucsc_color(color) for color in colors]
#         else:
#             colors = [
#                 matplotlib.colors.to_rgb(color) for color in re.split(r"\s|,|;", colors)
#             ]
#             return [color for i, color in zip(items, itertools.cycle(colors))]

#     bigwigs, statistics = infiles
#     sample_key = lambda d: d["track"].split(".")[0]

#     # Need to make an assembly hub if this is a custom genome
#     if not P.PARAMS["genome_custom"]:
#         hub, genomes_file, genome, trackdb = trackhub.default_hub(
#             hub_name=P.PARAMS["hub_name"],
#             short_label=P.PARAMS["hub_short"]
#             if not is_none(P.PARAMS["hub_short"])
#             else P.PARAMS["hub_name"],
#             long_label=P.PARAMS["hub_long"]
#             if not is_none(P.PARAMS["hub_long"])
#             else P.PARAMS["hub_name"],
#             email=P.PARAMS["hub_email"],
#             genome=P.PARAMS["genome_name"],
#         )

#         for sample, bigwigs_grouped in itertools.groupby(
#             sorted(bigwigs, key=sample_key), key=sample_key
#         ):

#             bigwigs_grouped = list(bigwigs_grouped)

#             # Create an overlay track
#             overlay = trackhub.AggregateTrack(
#                 aggregate="transparentOverlay",
#                 visibility="full",
#                 tracktype="bigWig",
#                 viewLimits="-2:2",
#                 maxHeightPixels="8:80:128",
#                 showSubtrackColorOnUi="on",
#                 name=f"{sample}",
#             )

#             for bigwig, color in zip(bigwigs_grouped, get_colors(bigwigs)):
#                 name = trackhub.helpers.sanitize(os.path.basename(bigwig))
#                 track = trackhub.Track(
#                     name=name,  # track names can't have any spaces or special chars.
#                     source=os.path.join(
#                         trackhub.helpers.data_dir(), os.path.basename(bigwig)
#                     ),  # filename to build this track from
#                     visibility="full",  # shows the full signal
#                     color=color,
#                     autoScale="on",  # allow the track to autoscale
#                     tracktype="bigWig",  # required when making a track
#                 )

#                 track_sub = trackhub.Track(
#                     name=f"{name}_subtrack",  # track names can't have any spaces or special chars.
#                     source=os.path.join(
#                         trackhub.helpers.data_dir(), os.path.basename(bigwig)
#                     ),  # filename to build this track from
#                     visibility="full",  # shows the full signal
#                     color=color,
#                     autoScale="on",  # allow the track to autoscale
#                     tracktype="bigWig",  # required when making a track
#                 )

#                 trackdb.add_tracks(track)
#                 overlay.add_subtrack(track_sub)

#         # Finalise trackhub
#         trackhub.upload.upload_hub(
#             hub=hub, host=P.PARAMS["hub_url"], remote_dir=P.PARAMS["hub_dir"]
#         )

#     else:
#         raise NotImplementedError("Custom genome not yet supported")

# hub = trackhub.Hub(P.PARAMS["hub_name"],
#                    short_label=P.PARAMS["hub_short"] if not is_none(P.PARAMS["hub_short"]) else P.PARAMS["hub_name"],
#                    long_label=P.PARAMS["hub_long"] if not is_none(P.PARAMS["hub_long"]) else P.PARAMS["hub_name"],
#                    email=P.PARAMS["hub_email"])


#         i
#         "longLabel": P.PARAMS["hub_long"]
#         if P.PARAMS["hub_long"]
#         else P.PARAMS["hub_name"],
#         "genomesFile": "genomes.txt",
#         "email": P.PARAMS["hub_email"],


# def get_track_data(fn):

#     track_name = " ".join(os.path.basename(fn).replace("_", " ").split(".")[:-1])

#     track_dict = {
#         "track": fn,
#         "bigDataUrl": f'{P.PARAMS["hub_url"].rstrip("/")}/{(os.path.join(P.PARAMS["ASSEMBLY_DIR"], fn)).lstrip("/")}',
#         "shortLabel": track_name,
#         "longLabel": track_name,
#         "type": f'{fn.split(".")[-1]}',
#     }

#     if P.PARAMS["hub_track_options"]:
#         try:
#             options = [op.strip() for op in P.PARAMS["hub_track_options"].split()]
#             options_dict = dict(zip(options[0::2], options[1::2]))
#             track_dict.update(options_dict)
#         except Exception as e:
#             print("Invalid custom track options")

#     return track_dict


# @mkdir(P.PARAMS["HUB_DIR"])
# @originate(os.path.join(P.PARAMS["HUB_DIR"], "hub.txt"))
# def generate_hub_metadata(outfile):

#     content = {
#         "hub": P.PARAMS["hub_name"],
#         "shortLabel": P.PARAMS["hub_short"]
#         if P.PARAMS["hub_short"]
#         else P.PARAMS["hub_name"],
#         "longLabel": P.PARAMS["hub_long"]
#         if P.PARAMS["hub_long"]
#         else P.PARAMS["hub_name"],
#         "genomesFile": "genomes.txt",
#         "email": P.PARAMS["hub_email"],
#         "descriptionUrl": "/".join(
#             [
#                 P.PARAMS["hub_url"].rstrip("/"),
#                 P.PARAMS["ASSEMBLY_DIR"].rstrip("/"),
#                 "visualise_statistics.html",
#             ]
#         ),
#     }

#     write_dict_to_file(outfile, content)


# @follows(generate_hub_metadata)
# @originate(os.path.join(P.PARAMS["HUB_DIR"], "genomes.txt"))
# def generate_assembly_metadata(outfile):

#     content = {
#         "genome": P.PARAMS["genome_name"],
#         "trackDb": os.path.join(P.PARAMS["genome_name"], "trackDb.txt"),
#     }

#     write_dict_to_file(outfile, content)


# @mkdir(P.PARAMS["ASSEMBLY_DIR"])
# @follows(generate_hub_metadata)
# @merge(make_bigwig, f'{P.PARAMS["ASSEMBLY_DIR"]}/trackDb.txt')
# def generate_trackdb_metadata(infiles, outfile):

#     # Generate all separate tracks
#     bigwig_tracks_all = [get_track_data(os.path.basename(fn)) for fn in infiles]

#     # Add colours to tracks
#     if not P.PARAMS["hub_colors"]:
#         colors = sns.color_palette("husl", len(bigwig_tracks_all))
#         for track, color in zip(bigwig_tracks_all, colors):
#             track["color"] = ",".join([str(c * 255) for c in color])
#     else:
#         for track, color in zip(
#             bigwig_tracks_all, itertools.cycle(P.PARAMS["hub_colors"].split(" "))
#         ):

#             track["color"] = ",".join(
#                 [str(c * 255) for c in matplotlib.colors.to_rgb(color)]
#             )

#     # Write track data separated
#     with open(outfile, "w") as w:
#         for track in bigwig_tracks_all:
#             for label, data in track.items():
#                 w.write(f"{label} {data}\n")
#             # Need to separate each track with a new line
#             w.write("\n")

#         # Group tracks by sample name and make separate combined tracks for each
#         sample_key = lambda d: d["track"].split(".")[0]
#         bigwig_tracks_grouped = {
#             sample: list(track)
#             for sample, track in itertools.groupby(
#                 sorted(bigwig_tracks_all, key=sample_key), key=sample_key
#             )
#         }

#         for sample, grouped_tracks in bigwig_tracks_grouped.items():

#             # Generate overlay track
#             combined_track_details = {
#                 "track": f"{sample}_combined",
#                 "container": "multiWig",
#                 "aggregate": "transparentOverlay",
#                 "showSubtrackColorOnUi": "on",
#                 "type": "bigWig 0 250",
#                 "shortLabel": f"{sample}_combined",
#                 "longLabel": f"{sample}_all_capture_probes_combined",
#                 "autoScale": "on",
#                 "windowingFunction": "maximum",
#             }

#             # Write overlay track
#             for label, data in combined_track_details.items():
#                 w.write(f"{label} {data}\n")
#             w.write("\n")

#             # Write sub-tracks
#             for track in grouped_tracks:
#                 track["track"] = track["track"].replace(".bigWig", "_subtrack.bigWig")
#                 for label, data in track.items():
#                     w.write(f"\t{label} {data}\n")
#                 w.write(f'\tparent {combined_track_details["track"]}\n')
#                 # Need to separate each track with a new line
#                 w.write("\n")


# @follows(generate_trackdb_metadata)
# @transform(
#     [make_bigwig, build_report], regex(".*/(.*)$"), P.PARAMS["ASSEMBLY_DIR"] + r"/\1"
# )
# def link_files(infile, outfile):
#     try:
#         infile_fp = os.path.abspath(infile)
#         os.symlink(infile_fp, outfile)
#     except Exception as e:
#         print(e)


# @originate("hub_url.txt")
# def write_hub_path(outfile):

#     with open(outfile, "w") as w:
#         w.write(
#             f'{P.PARAMS["hub_url"].rstrip("/")}/{os.path.join(P.PARAMS["HUB_DIR"], "hub.txt").lstrip("/")}\n'
#         )


def set_up_pipeline_params_dict():

    P.PARAMS["cluster_queue_manager"] = P.PARAMS.get("pipeline_cluster_queue_manager")
    P.PARAMS["conda_env"] = P.PARAMS.get(
        "conda_env", os.path.basename(os.environ["CONDA_PREFIX"])
    )


def set_up_chromsizes():

    if not is_none(P.PARAMS["genome_chrom_sizes"]):
        pass

    elif os.path.exists("chrom_sizes.txt.tmp"):
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"

    else:
        get_chromsizes_from_ucsc(P.PARAMS["genome_name"], "chrom_sizes.txt.tmp")
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"


@cli.command(
    context_settings=dict(
        ignore_unknown_options=True,
    )
)
@click.option("-h", "--help", default=False, is_flag=True)
@click.argument("mode", type=click.Choice(["make", "show", "clone", "touch"]))
# @click.option('--config', help='Path to pipeline configuration file', default='config.yml')
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def pipeline(mode, pipeline_options, help=False):

    sys.argv = ["pipeline", mode, *pipeline_options]

    if help:  # Don't want to waste time with set up if just want help with commands
        sys.argv.append(" --help")
    else:
        set_up_pipeline_params_dict()
        set_up_chromsizes()
        P.main()  # This runs the pipeline





# @follows(
#     mkdir("ccanalysis/bedgraph_compare"),
# )  # make_bedgraph)
# @collate(
#     "ccanalysis/bedgraphs/*.norm.bedgraph.gz",
#     regex(r".*/(?P<sample>.*)_(?:\d+)\.(?P<capture>.*)\.norm.bedgraph.gz"),
#     r"ccanalysis/bedgraph_compare/\1.\2.mean.bedgraph.gz",
# )
# def average_replicate_bedgraphs(infiles, outfile):

#     dframes = [
#         pd.read_csv(
#             fn,
#             sep="\t",
#             header=None,
#             chunksize=2e6,
#             names=["chrom", "start", "end", "score"],
#         )
#         for fn in infiles
#     ]

#     dframes_aggregated = []
#     for frames in zip(*dframes):
#         df_agg = (
#             pd.concat([d["score"] for d in frames], axis=1)
#             .assign(mean=lambda df: df.mean(axis=1), std=lambda df: df.std(axis=1))
#             .fillna(0)
#         )
#         df_agg = frames[0].loc[:, ["chrom", "start", "end"]].join(df_agg)
#         dframes_aggregated.append(df_agg)

#     df = pd.concat(dframes_aggregated)
#     # Only save the positive bins to save space
#     df_bdg_summary = df.query("mean > 0")[["mean", "std"]].to_csv(
#         outfile.replace(".mean.bedgraph.gz", ".tsv.gz")
#     )

#     # Output mean and union
#     df_bdg_mean = df[["chrom", "start", "end", "mean"]].to_csv(
#         outfile, sep="\t", index=None, header=None
#     )
#     df_bdg_union = df.loc[:, ~df.columns.isin(["mean", "std"])].to_csv(
#         outfile.replace(".mean.bedgraph.gz", ".union.bedgraph.gz")
#     )
