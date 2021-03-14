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
import os
import re
import sys
import pickle
from cgatcore import pipeline as P
from cgatcore.iotools import touch_file, zap_file

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
    originate
)
from ccanalyser.tools.statistics import (
    collate_slice_data,
    collate_read_data,
    collate_cis_trans_data,
    collate_histogram_data,
    extract_trimming_stats,
)

from ccanalyser.utils import is_on, is_none, make_group_track


##############################
#   Set-up global parameters #
##############################

P.get_parameters("config.yml")

##############################
#  Pipeline stages           #
##############################


def set_up_pipeline_params_dict():

    P.PARAMS["cluster_queue_manager"] = P.PARAMS.get("pipeline_cluster_queue_manager")
    P.PARAMS["conda_env"] = P.PARAMS.get(
        "conda_env", os.path.basename(os.environ["CONDA_PREFIX"])
    )

    # Sanitise hub name
    P.PARAMS["hub_name"] = re.sub(r"[,\s+\t;:]", "_", P.PARAMS["hub_name"])


def set_up_chromsizes():

    if not is_none(P.PARAMS["genome_chrom_sizes"]):
        pass

    elif os.path.exists("chrom_sizes.txt.tmp"):
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"

    else:
        from pybedtools.helpers import get_chromsizes_from_ucsc
        get_chromsizes_from_ucsc(P.PARAMS["genome_name"], "chrom_sizes.txt.tmp")
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"

@follows(mkdir("pre_ccanalysis"), mkdir("pre_ccanalysis/fastqc"))
@transform("*.fastq*", regex(r"(.*).fastq.*"), r"pre_ccanalysis/fastqc/\1_fastqc.zip")
def qc_reads(infile, outfile):
    
    """Runs fastqc on the input files"""
    
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

    tmp = outfile.replace(".gz", "")
    statement = """ccanalyser genome digest
                   %(infile)s
                   -l %(tmp)s.log
                   -o %(tmp)s
                   -r %(analysis_restriction_enzyme)s
                   --sort
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
    """Splits the input fastq files into chunks for parallel processing"""

    infiles = " ".join(infiles)
    output_prefix = outfile.replace(".log", "")

    statement = """ccanalyser fastq split
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
    extras=[r"\1", r"\2"],
)
def fastq_duplicates_parse(infiles, outfile, sample_name, part_no):

    """Parses fastq files into json format for sequence deduplication"""

    fq1, fq2 = [os.path.abspath(fn) for fn in infiles]

    statement = """ccanalyser fastq deduplicate parse
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
    fastq_duplicates_parse,
    regex(r"pre_ccanalysis/deduplicated/deduplicated_ids/(.*)_\d*.json.gz"),
    r"pre_ccanalysis/deduplicated/deduplicated_ids/\1.json.gz",
)
def fastq_duplicates_identify(infiles, outfile):

    """Identifies duplicate sequences from parsed fastq files (json format)"""

    infiles_str = " ".join(infiles)
    statement = """ccanalyser fastq deduplicate identify
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
    fastq_duplicates_parse,
    fastq_duplicates_identify,
    mkdir("statistics/deduplication/data/"),
)
@collate(
    "pre_ccanalysis/split/*.fastq*",
    regex(r".*/(.*_part\d+)_[12].fastq(?:.gz)?"),
    r"pre_ccanalysis/deduplicated/\1_1.fastq",
)
def fastq_duplicates_remove(infiles, outfile):

    """Removes duplicate read fragments identified by fastq_duplicates_identify
       from input fastq files"""

    fq1, fq2 = infiles
    sample = re.match(r".*/(.*)(_part\d+)_[12].fastq(?:.gz)?", fq1)
    sample_name = sample.group(1)
    sample_part = sample.group(2)
    dd_ids = f"pre_ccanalysis/deduplicated/deduplicated_ids/{sample_name}.json.gz"
    stats_prefix = f"statistics/deduplication/data/{sample_name}_{sample_part}"
    output_prefix = outfile.replace("_1.fastq", "")

    statement = """ccanalyser fastq deduplicate remove
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


@follows(fastq_duplicates_remove)
@merge(
    "statistics/deduplication/data/*.csv",
    "statistics/deduplication/deduplication.reads.csv",
)
def collate_deduplication_stats(infiles, outfile):

    """Combines deduplication statistics from partitions."""

    stats_prefix = outfile.replace(".reads.csv", "")

    df_stats = collate_read_data(infiles)

    df_stats_read = df_stats.query('stat_type != "reads_removed"')

    df_stats.to_csv(f"{stats_prefix}.summary.csv", index=False)
    df_stats_read.to_csv(
        outfile, index=False
    )  # Modified to enable more streamlined summary at final stage


@follows(
    mkdir("pre_ccanalysis/trimmed"),
    fastq_duplicates_remove,
    mkdir("statistics/trimming/data/"),
)
@collate(
    "pre_ccanalysis/deduplicated/*.fastq*",
    regex(r"pre_ccanalysis/deduplicated/(.*)_[12].fastq(?:.gz)?"),
    r"pre_ccanalysis/trimmed/\1_1_val_1.fq.gz",
)
def trim_reads(infiles, outfile):
    
    """Trim adaptor sequences from fastq files using Trim-galore"""

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

    """Extracts and collates adapter trimming statistics"""

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
    sn = fn.split("_part")[0]
    statement = """ccanalyser fastq digest
                   %(infile)s
                   -m flashed
                   -r %(analysis_restriction_enzyme)s
                   -o %(outfile)s
                   --stats_prefix %(outfile)s.stats
                   --minimum_slice_length 18
                   --compression_level %(pipeline_advanced_compression)s
                   -p 1
                   --stats_prefix statistics/digestion/data/%(fn)s.flashed
                   --sample_name %(sn)s
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
    sn = fn.split("_part")[0]

    statement = """ccanalyser fastq digest
                   %(fq1)s
                   %(fq2)s
                   -m pe
                   -r %(analysis_restriction_enzyme)s
                   -o %(outfile)s
                   --minimum_slice_length 18
                   -p 1
                   --compression_level %(pipeline_advanced_compression)s
                   --stats_prefix statistics/digestion/data/%(fn)s.pe
                   --sample_name %(sn)s
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

    """Combines digestion statistics"""

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

    '''Indexes all bam files (both partitioned and merged)'''

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


@originate("ccanalysis/annotations/exclude.bed")
def build_exclusion_bed(outfile):
    
    """Generates exclusion window around each capture site"""

    statement = """bedtools slop
                    -i %(analysis_capture_oligos)s -g %(genome_chrom_sizes)s -b %(analysis_reporter_exclusion_zone)s
                    | bedtools subtract -a - -b %(analysis_capture_oligos)s
                    | sort -k1,1 -k2,2n 
                    > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@originate("ccanalysis/annotations/capture.bed")
def sort_capture_oligos(outfile):

    """Sorts the capture oligos for bedtools intersect with --sorted option"""

    statement = """cat %(analysis_capture_oligos)s | sort -k1,1 -k2,2n > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(P.PARAMS.get("analysis_blacklist"))
@originate("ccanalysis/annotations/blacklist.bed")
def sort_blacklist(outfile):

    """Sorts the capture oligos for bedtools intersect with --sorted option"""

    if os.path.exists(P.PARAMS['analysis_blacklist']):
        statement = """cat %(analysis_blacklist)s | sort -k1,1 -k2,2n > %(outfile)s"""
    else:
        statement = 'touch %(outfile)s'
    
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
                "fn": "ccanalysis/annotations/capture.bed",
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
                "fn": "ccanalysis/annotations/capture.bed",
                "action": "count",
                "fraction": 0.9,
            },
            {
                "name": "blacklist",
                "fn": "ccanalysis/annotations/blacklist.bed",
                "action": "count",
                "fraction": 1e-9,
            },
        ]
    ),
    r"ccanalysis/annotations/\1.annotations.tsv",
)
def annotate_slices(infile, outfile):

    '''Annotates mapped read slices.
     
       Slices are annotated with:
       * capture name 
       * capture count
       * exclusion name
       * exclusion count
       * blacklist count
       * restriction fragment number
       '''

    slices = infile[0]
    flags = {"name": "-n", "fn": "-b", "action": "-a", "fraction": "-f"}

    cmd_args = []
    for args in infile[1]:
        for arg_name, arg in args.items():
            cmd_args.append(f'{flags.get(arg_name)} {arg if arg else "-"}')

    cmd_args = " ".join(cmd_args)
    statement = """ccanalyser alignments annotate
                    %(slices)s
                    %(cmd_args)s
                    -o %(outfile)s
                    --invalid_bed_action ignore
                    -p 4
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
    """Runs the pipeline until just prior to identification of reporters"""
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
def identify_reporters(infiles, outfile):
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

    statement = """ccanalyser
                   reporters
                   identify
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


@follows(mkdir("ccanalysis/reporters/combined"), identify_reporters)
@collate(
    "ccanalysis/reporters/unfiltered/*.tsv",
    regex(
        r".*/(?P<sample>.*)_part\d+.(flashed|pe).(?P<capture>.*).(slices|fragments).tsv"
    ),
    r"ccanalysis/reporters/combined/\1.\2.\3.\4.tsv",
    extras=[r"\1", r"\2", r"\3", r"\4"],
)
def collate_reporters(infiles, outfile, *grouping_args):

    """Concatenates identified reporters """

    # Need to concat tsv files but remove headers, the sed command performs the removal
    statement = (
        f"cat {' '.join(infiles)} | sed -e '1n' -e '/.*parent_read.*/d' > {outfile}"
    )

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zero un-aggregated reporters
    for fn in infiles:
        zap_file(fn)


@follows(identify_reporters, mkdir("ccanalysis/reporters/deduplicated"))
@transform(
    collate_reporters,
    regex(r".*/(?P<sample>.*).(flashed|pe).(?P<capture>.*).fragments.tsv"),
    r"ccanalysis/reporters/deduplicated/\1.\2.\3.json.gz",
    extras=[r"\2"],
)
def deduplicate_fragments(infile, outfile, read_type):

    '''Identifies duplicate fragments by removing fragments with the same coordinates
       and slice order'''

    statement = """ccanalyser alignments deduplicate
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
    extras=[r"\1", r"\2", r"\3"],
)
def deduplicate_reporters(infile, outfile, sample_name, read_type, capture_oligo):

    '''Removes reporters with duplicate coordinates'''

    slices, duplicated_ids = infile
    stats_prefix = (
        f"statistics/ccanalysis/data/{sample_name}_{read_type}_{capture_oligo}"
    )

    statement = """ccanalyser alignments deduplicate
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
    deduplicate_reporters,
    regex(r".*/(?P<sample>.*).(?:flashed|pe).(?P<capture>.*).slices.tsv"),
    r"ccanalysis/reporters/\1.\2.tsv.gz",
    extras=[r"\1", r"\2"],
)
def collate_reporters_final(infiles, outfile, *grouping_args):

    '''Final collation of reporters by sample and capture probe'''

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

    ''''Combination of all reporter identification and filtering statistics'''

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


@follows(deduplicate_reporters, collate_ccanalyser_stats)
def post_ccanalysis():
    '''Reporters have been identified, deduplicated and collated by sample/capture probe'''



@follows(mkdir("ccanalysis/interactions/counts"))
@transform(
    collate_reporters_final,
    regex(r"ccanalysis/reporters/(.*)\.(.*).tsv.gz"),
    r"ccanalysis/interactions/counts/\1.\2.tsv.gz",
)
def count_interactions(infile, outfile):

    '''Counts the number of interactions identified between reporter restriction fragments'''

    statement = [
        "ccanalyser reporters  count",
        "%(infile)s",
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

@follows(mkdir('ccanalysis/interactions/fragments'))
@transform(
    count_interactions,
    regex(r"ccanalysis/interactions/counts/(.*)\.(.*)\.tsv.gz"),
    add_inputs(digest_genome),
    r"ccanalysis/interactions/fragments/\1.\2.fragments.hdf5",
    extras=[r"\1", r"\2"],
)
def store_interactions_at_fragment_level(infile, outfile, sample_name, capture_name):

    '''Stores restriction fragment interaction counts in cooler format'''

    counts, rf_map = infile
    output_prefix = outfile.replace(f".{capture_name}.fragments", "")

    statement = """ccanalyser reporters  store
                   fragments
                   %(counts)s
                   -f %(rf_map)s
                   -g %(genome_name)s
                   -n %(capture_name)s
                   -o %(output_prefix)s
                   -c %(analysis_capture_oligos)s
                   --suffix fragments
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@follows(digest_genome, count_interactions)
@originate(r'ccanalysis/interactions/binners.pkl')
def generate_bin_conversion_tables(outfile):
    """Converts restriction fragments to genomic bins."
    """

    from ccanalyser.tools.storage import GenomicBinner    

    frags = pd.read_csv('pre_ccanalysis/restriction_enzyme_map/genome.digest.bed.gz', 
                        sep='\t', 
                        names=['chrom', 'start', 'end', 'name'])
    
    binner_dict = dict()
    for bs in re.split(r"[,;]\s*|\s+", str(P.PARAMS["plot_bin_size"])):
        gb = GenomicBinner(chromsizes=P.PARAMS['genome_chrom_sizes'], fragments=frags, binsize=int(bs))
        gb.bin_conversion_table # Property that is cached so need to call it to make sure it is present.
        binner_dict[int(bs)] = gb
   
    with open('ccanalysis/interactions/binners.pkl', 'wb') as w:
        pickle.dump(binner_dict, w)



@active_if(not is_none(P.PARAMS.get("plot_bin_size")))
@follows(generate_bin_conversion_tables, mkdir('ccanalysis/interactions/binned/'))
@transform(
    store_interactions_at_fragment_level,
    regex(r"ccanalysis/interactions/fragments/(.*)\.(.*)\.fragments\.hdf5"),
    add_inputs(generate_bin_conversion_tables),
    r"ccanalysis/interactions/binned/\1.\2.log",
    extras=[r"\2"],
)
def store_interactions_binned(infile, outfile, capture_name):

    '''
    Converts a cooler file of restriction fragments to even genomic bins.
    '''

    infile, conversion_tables = infile

    bin_options = " -b " + " -b ".join(
        re.split(r"[,;]\s*|\s+", str(P.PARAMS["plot_bin_size"]))
    )
    output_prefix = outfile.replace(f".{capture_name}.log", "")

    statement = f"""ccanalyser reporters  store
                  bins
                  %(infile)s
                  -o %(output_prefix)s
                  %(bin_options)s
                  --conversion_tables %(conversion_tables)s
                  --normalise
                  -p 4
                  > %(outfile)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=4,
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(store_interactions_at_fragment_level, store_interactions_binned)
@collate(
    ["ccanalysis/interactions/fragments/*.hdf5", "ccanalysis/interactions/binned/*.hdf5"],
    regex(r".*/(.*)\.(.*)\.(?:fragments|\d+)\.hdf5"),
    r"ccanalysis/interactions/\1.hdf5",
    extras=[r"\1"],
)
def merge_interactions(infiles, outfile, sample_name):

    '''Combines cooler files together'''

    infiles_str = " ".join(infiles)
    statement = f"""ccanalyser reporters  store
                     merge
                     %(infiles_str)s
                     -o %(outfile)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    for fn in infiles:
        zap_file(fn)
    
    zap_file('ccanalysis/interactions/binners.pkl')
    

@follows(mkdir("ccanalysis/bedgraphs"))
@transform(
    merge_interactions,
    regex(r".*/(.*).hdf5"),
    r"ccanalysis/bedgraphs/\1.raw.log",
    extras=[r"\1"],
)
def make_bedgraph_raw(infile, outfiles, sample_name):
    """Extract reporters in bedgraph format from stored interactions"""

    output_prefix = f"ccanalysis/bedgraphs/{sample_name}.raw"

    statement = """ccanalyser reporters  bedgraph
                   %(infile)s
                   -o %(output_prefix)s
                   > %(output_prefix)s.log"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(
    merge_interactions,
    regex(r".*/(.*).hdf5"),
    r"ccanalysis/bedgraphs/\1.normalised.log",
    extras=[r"\1"],
)
def make_bedgraph_normalised(infile, outfiles, sample_name):
    """Extract reporters in bedgraph format from stored interactions.
       Normalises the counts by the number of cis interactions identified."""

    output_prefix = f"ccanalysis/bedgraphs/{sample_name}.normalised"

    statement = """ccanalyser reporters  bedgraph
                   %(infile)s
                   -o %(output_prefix)s
                   --normalise
                   > %(output_prefix)s.log
                   """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


# @transform(
#     merge_interactions,
#     regex(r".*/(.*).hdf5"),
#     r"ccanalysis/bedgraphs/\1.windowed.log",
#     extras=[r"\1"],
# )
# def make_bedgraph_windowed(infile, outfiles, sample_name):
#     """Extract reporters in bedgraph format from stored interactions"""

#     output_prefix = f"ccanalysis/bedgraphs/{sample_name}.windowed"

#     statement = """ccanalyser reporters  bedgraph
#                    %(infile)s
#                    -o %(output_prefix)s
#                    --binsize 5000
#                    --normalise
#                    > %(output_prefix)s.log
#                    """
#     P.run(
#         statement,
#         job_queue=P.PARAMS["pipeline_cluster_queue"],
#         job_condaenv=P.PARAMS["conda_env"],
#     )


@follows(
    mkdir("ccanalysis/bigwigs"),
    make_bedgraph_raw,
    make_bedgraph_normalised,
)
@transform(
    "ccanalysis/bedgraphs/*",
    regex(r"ccanalysis/bedgraphs/(.*).bedgraph"),
    r"ccanalysis/bigwigs/\1.bigWig",
)
def make_bigwig(infile, outfile):
    """Uses UCSC tools bedGraphToBigWig to generate bigWigs for each bedgraph"""

    tmp = f"{outfile}.tmp"
    statement = """  cat %(infile)s
                   | sort -k1,1 -k2,2n > %(tmp)s
                   && bedGraphToBigWig %(tmp)s %(genome_chrom_sizes)s %(outfile)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("capture_compare/bedgraphs_union"))
@collate(
    "ccanalysis/bedgraphs/*.bedgraph",
    regex(r".*/(?:.*)\.(raw|normalised|windowed)\.(.*).bedgraph"),
    r"capture_compare/bedgraphs_union/\2.\1.tsv",
    extras=[r"\1", r"\2"],
)
def make_union_bedgraph(infiles, outfile, normalisation_type, capture_name):

    '''Collates bedgraphs by capture probe into a single file for comparison.'''

    infiles_str = " ".join(infiles)
    sample_names = " ".join(
        [
            re.match(r".*/(.*)\.(.*)\.(?:.*).bedgraph(?:.gz)?", fn).group(1)
            for fn in infiles
        ]
    )

    statement = """bedtools 
                  unionbedg 
                  -i %(infiles_str)s 
                  -header 
                  -names %(sample_names)s
                  > %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

@merge(
    [collate_deduplication_stats, collate_digestion_stats, collate_ccanalyser_stats],
    "statistics/run_statistics.csv",
)
def merge_run_statistics(infiles, outfile):

    '''Generates a summary statistics file for the pipeline run.
       Summarised at the read count level.'''

    df = pd.concat([pd.read_csv(fn) for fn in infiles])

    df.sort_values(
        ["sample", "read_type", "stat"], ascending=[True, True, False]
    ).to_csv(outfile)


@merge(
    [merge_run_statistics],
    "statistics/visualise_statistics.html",
)
def build_report(infile, outfile):
    """Run jupyter notebook for reporting and plotting pipeline statistics"""

    path_script = __file__
    path_script_dir = os.path.dirname(path_script)
    path_nb_dir = os.path.dirname(path_script_dir)

    statement = """rm statistics/visualise_statistics* -f &&
                   papermill
                   %(path_nb_dir)s/visualise_capture-c_stats.ipynb
                   statistics/visualise_statistics.ipynb
                   -p directory $(pwd)/statistics/ &&
                   jupyter nbconvert
                   --no-input
                   --to html
                   statistics/visualise_statistics.ipynb
                   statistics/visualise_statistics.html
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(is_on(P.PARAMS.get("hub_create")))
@collate(
    make_bigwig,
    regex(r".*/(?:.*)\.normalised\.(?:.*).bigWig"),
    os.path.join(P.PARAMS.get("hub_dir", ''), P.PARAMS.get("hub_name", '') + ".hub.txt"),
    extras=[build_report],
)
def make_ucsc_hub(infiles, outfile, statistics):
    '''Creates a ucsc hub from the pipeline output'''

    import trackhub

    bigwigs = infiles
    key_sample = lambda b: os.path.basename(b).split(".")[0]
    key_capture = lambda b: b.split(".")[-2]
    

    # Need to make an assembly hub if this is a custom genome
    if not P.PARAMS["genome_custom"]:
        hub, genomes_file, genome, trackdb = trackhub.default_hub(
            hub_name=P.PARAMS["hub_name"],
            short_label=P.PARAMS.get("hub_short"),
            long_label=P.PARAMS.get("hub_long"),
            email=P.PARAMS["hub_email"],
            genome=P.PARAMS["genome_name"],
        )
    
        for key in [key_sample, key_capture]:

            trackdb.add_tracks(make_group_track(bigwigs, key, overlay=True).values())

        
        if is_on(
            P.PARAMS.get("hub_upload")
        ):  # If the hub need to be uploaded to a server
            trackhub.upload.upload_hub(
                hub=hub, host=P.PARAMS["hub_url"], remote_dir=P.PARAMS["hub_dir"]
            )
        else:
            trackhub.upload.stage_hub(hub=hub, staging=P.PARAMS["hub_dir"])

    else:
        raise NotImplementedError("Custom genome not yet supported")


@follows(make_ucsc_hub)
@originate("hub_url.txt")
def write_hub_path(outfile):
    '''Convinence task to write hub url to use for adding custom hub to UCSC genome browser'''

    with open(outfile, "w") as w:
        url = P.PARAMS["hub_url"].rstrip("/")
        name_dir = P.PARAMS["hub_dir"].strip("/")
        name_hubtxt = P.PARAMS["hub_name"] + ".hub.txt"

        path_hubtxt = f"{url}/{name_dir}/{name_hubtxt}"

        w.write(path_hubtxt)


@follows(mkdir('capture_compare/differential'))
@transform('capture_compare/bedgraphs_union/*.tsv',
           regex(r'.*/(.*)\.raw\.tsv'),
           r'capture_compare/differential/\1.log',
           extras=[r'\1'])
def identify_differential_interactions(infile, outfile, capture_name):

    if len(infile) >= 4:

        output_prefix = outfile.replace('.log', '')

        statement = '''ccanalyser
                    interactions differential
                    %(infile)s
                    -n %(capture_name)s
                    -c %(analysis_capture_oligos)s
                    -o %(output_prefix)s
                    > %(outfile)s
                    '''
        
        P.run(
            statement,
            job_queue=P.PARAMS["pipeline_cluster_queue"],
            job_condaenv=P.PARAMS["conda_env"],
        )
    
    else:
        print('Not enough replicates for differential testing')

@follows(merge_interactions, mkdir("ccanalysis/heatmaps/"))
@transform(
    "ccanalysis/interactions/*.hdf5",
    regex(r"ccanalysis/interactions/(.*).hdf5"),
    r"ccanalysis/heatmaps/\1.log",
)
def plot_interactions(infile, outfile):
    """Plots a heatmap over a specified region"""

    #TODO: Plot multiple resolutions

    if not is_none(P.PARAMS.get("plot_normalisation")):
        norm = P.PARAMS['plot_normalisation']
    else:
        norm_default = {'capture': 'n_interactions',
                        'tri': "n_rf_n_interactions",
                        'tiled': "ice",}
        norm = norm_default[P.PARAMS['analysis_method']]
    
    
    output_prefix = outfile.replace('.log', '')

    resolutions = ' -r '.join(
        re.split(r"[,;]\s*|\s+", str(P.PARAMS["plot_bin_size"]))
    )

    statement = f"""ccanalyser reporters  plot
                    %(infile)s
                    -r %(resolutions)s
                    -c %(plot_coordinates)s
                    --normalisation %(norm)s
                    --cmap {P.PARAMS.get('plot_cmap', 'jet')}
                    --vmin {P.PARAMS.get('plot_min', '0')}
                    --vmax {P.PARAMS.get('plot_vmax', '1')}
                    -o %(output_prefix)s
                    > %(outfile)s 2>&1"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
        )

@merge([make_ucsc_hub, plot_interactions], 'pipeline_complete.txt')
def full(infiles, outfile):
    touch_file(outfile)


if __name__ == '__main__':
    set_up_chromsizes()
    set_up_pipeline_params_dict()
    P.main(sys.argv)


