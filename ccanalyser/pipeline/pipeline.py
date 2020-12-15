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

import itertools
import os
import re
import sys

import matplotlib.colors
import seaborn as sns
import trackhub
from ccanalyser.ccanalyser_cli import PACKAGE_DIR
from ccanalyser.utils.helpers import is_none, is_off, is_on
from cgatcore import pipeline as P
from pybedtools.helpers import get_chromsizes_from_ucsc
from ruffus import (active_if, add_inputs, collate, follows, merge, mkdir,
                    originate, regex, split, suffix, transform)

# Global vars
# TODO: Issue with selecting queue manager -- defaults to sun grid (cgat-core issue)
## kwargs.get("cluster_queue_manager") is the offending option

# Read in parameter file
params = "config.yml"
P.get_parameters(params)

# Sort pipeline global params
P.PARAMS["cluster_queue_manager"] = P.PARAMS["pipeline_cluster_queue_manager"]
P.PARAMS["conda_env"] = P.PARAMS.get(
    "conda_env", os.path.basename(os.environ["CONDA_PREFIX"])
)
# P.PARAMS['HUB_DIR'] = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
# P.PARAMS['ASSEMBLY_DIR'] = os.path.join(P.PARAMS['HUB_DIR'], P.PARAMS['genome_name'])

# Check if chromsizes are provided, if not download them
if is_none(P.PARAMS["genome_chrom_sizes"]):
    get_chromsizes_from_ucsc(P.PARAMS["genome_name"], "chrom_sizes.txt.tmp")
    P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"


@follows(mkdir("fastq"), mkdir("fastq/fastqc"))
@transform("*.fastq*", regex(r"(.*).fastq.*"), r"fastq/fastqc/\1_fastqc.zip")
def qc_reads(infile, outfile):
    """Quality control of raw sequencing reads"""
    outdir = os.path.dirname(outfile)
    statement = """fastqc
                   -q
                   -t %(pipeline_n_cores)s
                   --nogroup %(infile)s
                   --outdir %(outdir)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("run_statistics"))
@merge(qc_reads, "run_statistics/fastqc_report.html")
def multiqc_reads(infile, outfile):
    """Collate fastqc reports into single report using multiqc"""

    bn = os.path.basename(outfile)
    dn = os.path.dirname(outfile)

    statement = """rm -f %(outfile)s &&
                   export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc
                   fastq/fastqc/
                   -o %(dn)s
                   -n %(bn)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("ccanalysis/restriction_enzyme_map/"))
@transform(
    P.PARAMS["genome_fasta"],
    regex(r".*/(.*).fa.*"),
    r"ccanalysis/restriction_enzyme_map/genome.digest.bed.gz",
)
def digest_genome(infile, outfile):
    """Digest genome using restriction enzyme and output fragments in bed file"""

    tmp = outfile.replace(".gz", "")
    restriction_enzyme = P.PARAMS["analysis_restriction_enzyme"]

    statement = """ccanalyser utils digest_genome
                   -i %(infile)s
                   -l %(tmp)s.log
                   -o %(tmp)s
                   -r %(restriction_enzyme)s
                   && gzip %(tmp)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("fastq/split"))
@collate(
    "*.fastq.gz", regex(r"(.*)_R*[12].fastq.*"), r"fastq/split/\1_part0_1.fastq.gz",
)
def split_fastq(infiles, outfile):
    """Splits the combined (flashed) fastq files into chunks for parallel processing"""
    
    infiles = ' '.join(infiles)
    output_prefix = outfile.replace("_part0_1.fastq.gz", "")

    statement = """ccanalyser utils split_fastq
                  %(infiles)s
                  -o %(output_prefix)s
                  -n %(split_n_reads)s
                  -c %(pipeline_advanced_compression)s """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(is_on(P.PARAMS["deduplication_pre-dedup"]))
@follows(mkdir("fastq/deduplicated"), mkdir("fastq/deduplicated/deduplicated_ids"), split_fastq)
@collate(
    'fastq/split/*.fastq.gz',
    regex(r"fastq/split/(.*)_part(\d+)_[12].fastq.gz"),
    r"fastq/deduplicated/deduplicated_ids/\1_\2.pkl",
)
def find_duplicate_reads(infiles, outfile):

    """Checks for duplicate read1/read2 pairs in a pair of fastq files
    any duplicates are discarded"""

    fq1, fq2 = [os.path.abspath(fn) for fn in infiles]

    statement = """ccanalyser utils deduplicate_fastq find_duplicates
                            %(fq1)s %(fq2)s
                            -d %(outfile)s
                            -c %(pipeline_advanced_compression)s
                            """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="8G",
        job_condaenv=P.PARAMS["conda_env"],
    )

@collate(find_duplicate_reads,
         regex(r'fastq/deduplicated/deduplicated_ids/(.*)_\d*.pkl'),
         r'fastq/deduplicated/deduplicated_ids/\1.pkl')
def merge_read_ids(infiles, outfile):
    infiles = ' '.join(infiles)
    statement = """ccanalyser utils deduplicate_fastq merge_ids
                            %(infiles)s
                            -o %(outfile)s
                            """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="8G",
        job_condaenv=P.PARAMS["conda_env"],
    )

@follows(find_duplicate_reads, merge_read_ids)
@collate(
    'fastq/split/*.fastq.gz',
    regex(r"fastq/split/(.*\d)_[12].fastq.gz"),
    add_inputs(merge_read_ids),
    r"fastq/deduplicated/\1_1.fastq.gz",
    )
def remove_duplicate_reads(infiles, outfile):

    """Checks for duplicate read1/read2 pairs in a pair of fastq files
    any duplicates are discarded"""

    fq_files, dd_ids = zip(*infiles)
    fq1, fq2 = fq_files

    # Issue adding the correct merged id file so will select the correct one here
    fq_original_prefix = outfile.split('/')[-1].split('_part')[0]
    dd_ids = [ids for ids in dd_ids if f'{fq_original_prefix}.pkl' in ids][0]
    out1, out2 = outfile, outfile.replace("_1.fastq.gz", "_2.fastq.gz")


    statement = """ccanalyser utils deduplicate_fastq remove_duplicates
                            %(fq1)s %(fq2)s
                            -d %(dd_ids)s
                            -o %(out1)s %(out2)s
                            -c %(pipeline_advanced_compression)s
                            """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="8G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("fastq/trimmed"), remove_duplicate_reads)
@collate(
    r"fastq/deduplicated/*.fastq.gz",
    regex(r"fastq/deduplicated/(.*)_[12].fastq.gz"),
    r"fastq/trimmed/\1_1_val_1.fq.gz",
)
def trim_reads(infiles, outfile):
    """Trim adaptor sequences using Trim-galore"""

    fastq1, fastq2 = infiles
    outdir = os.path.dirname(outfile)
    statement = """trim_galore
                   --cores %(pipeline_n_cores)s
                   --paired %(trim_options)s
                   -o %(outdir)s
                   %(fastq1)s
                   %(fastq2)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(trim_reads, mkdir("fastq/flashed"))
@collate(
    "fastq/trimmed/*.fq.gz",
    regex(r"fastq/trimmed/(.*)_[12]_.*.fq.gz"),
    r"fastq/flashed/\1.extendedFrags.fastq.gz",
)
def combine_reads(infiles, outfile):
    """Combine overlapping paired-end reads using flash"""
    fastq1, fastq2 = infiles
    output_prefix = outfile.replace(".extendedFrags.fastq.gz", "")
    statement = """flash
                   -z
                   -t %(pipeline_n_cores)s
                   -o %(output_prefix)s
                    %(fastq1)s
                    %(fastq2)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("fastq/digested"), combine_reads)
@transform(
    'fastq/flashed/*.fastq.gz',
    regex(r"fastq/flashed/(.*).extendedFrags.fastq.gz"),
    r"fastq/digested/\1.flashed.fastq.gz",
)
def digest_flashed_reads(infile, outfile):
    """In silico restriction enzyme digest of combined (flashed) read pairs"""

    statement = """ccanalyser utils digest_fastq
                   flashed
                   -o %(outfile)s
                   --stats_file %(outfile)s.stats
                   -r %(analysis_restriction_enzyme)s
                   -m 18
                   -c %(pipeline_advanced_compression)s
                   -p 1
                   -i %(infile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=4,
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(combine_reads)
@collate(
    'fastq/flashed/*.fastq.gz',
    regex(r"fastq/flashed/(.*).notCombined_[12].fastq.gz"),
    r"fastq/digested/\1.pe.fastq.gz",
)
def digest_pe_reads(infiles, outfile):
    """In silico restriction enzyme digest of non-combined (non-flashed) read pairs"""

    fq1, fq2 = infiles

    statement = """ccanalyser utils digest_fastq
                   unflashed
                   --stats_file %(outfile)s.stats
                   -r %(analysis_restriction_enzyme)s
                   -o %(outfile)s
                   -c %(pipeline_advanced_compression)s
                   -m 18
                   -p 1
                   -1 %(fq1)s
                   -2 %(fq2)s
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=4,
        job_condaenv=P.PARAMS["conda_env"],
    )

@follows(multiqc_reads,
         digest_flashed_reads,
         digest_pe_reads)
def fastq_preprocessing():
    pass


@follows(mkdir("aligned"))
@transform(
    [digest_flashed_reads, digest_pe_reads],
    regex(r"fastq/digested/(.*).fastq.gz"),
    r"aligned/\1.bam",
)
def align_reads(infile, outfile):
    """ Aligns digested fq files using bowtie2"""
    aligner = P.PARAMS["align_aligner"]
    index_flag = (
        P.PARAMS["align_index_flag"]
        if not is_none(P.PARAMS["align_index_flag"])
        else ""
    )
    options = (
        P.PARAMS["align_options"] if not is_none(P.PARAMS["align_options"]) else ""
    )

    statement = """%(aligner)s %(options)s %(index_flag)s %(genome_aligner_index)s %(infile)s
                    | samtools view -b -S > %(outfile)s
                    && samtools sort %(outfile)s -o %(outfile)s.sorted.bam -m 2G -@ %(pipeline_n_cores)s
                    && mv -f %(outfile)s.sorted.bam %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory="4G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@collate(align_reads, regex(r"aligned/(.*)_part\d+.*.bam"), r"aligned/\1.bam")
def merge_bam_files(infiles, outfile):
    """Combines bam files (by flashed/non-flashed status and sample)"""
    fnames = " ".join(infiles)

    statement = """samtools merge %(outfile)s %(fnames)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(merge_bam_files, regex("(.*).bam"), r"\1.bam.bai")
def index_bam(infile, outfile):
    statement = """samtools index %(infile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="1G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("aligned/mapping_statistics"), index_bam)
@transform(
    merge_bam_files,
    regex(r"aligned/(.*).bam"),
    r"aligned/mapping_statistics/\1.picard.metrics",
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


@merge(mapping_qc, "run_statistics/mapping_report.html")
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
    align_reads, regex(r"aligned/(.*).bam"), r"ccanalysis/annotations/\1.bam.bed.gz"
)
def bam_to_bed(infile, outfile):
    """Converts bam files to bed for faster intersection"""
    tmp = outfile.replace(".gz", "")
    statement = """bedtools bamtobed -i %(infile)s | gzip > %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(
    P.PARAMS["analysis_capture_oligos"],
    regex(P.PARAMS["analysis_capture_oligos"]),
    r"ccanalysis/annotations/exclude.bed.gz",
)
def build_exclusion_bed(infile, outfile):
    """Generates exclusion window around each capture site"""

    statement = """bedtools slop
                    -i %(infile)s -g %(genome_chrom_sizes)s -b %(analysis_reporter_exclusion_zone)s
                    | bedtools subtract -a - -b %(analysis_capture_oligos)s
                    | gzip > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(digest_genome, build_exclusion_bed)
@transform(
    bam_to_bed,
    regex(r"ccanalysis/annotations/(.*).bam.bed.gz"),
    add_inputs(
        [
            {
                "name": "restriction_fragment",
                "fn": "ccanalysis/restriction_enzyme_map/genome.digest.bed.gz",
                "action": "get",
                "overlap_fraction": 0.2,
            },
            {
                "name": "capture",
                "fn": P.PARAMS["analysis_capture_oligos"],
                "action": "get",
                "overlap_fraction": 0.9,
            },
            {
                "name": "exclusion",
                "fn": "ccanalysis/annotations/exclude.bed.gz",
                "action": "get",
                "overlap_fraction": 1e-9,
            },
            {
                "name": "exclusion_count",
                "fn": "ccanalysis/annotations/exclude.bed.gz",
                "action": "count",
                "overlap_fraction": 1e-9,
            },
            {
                "name": "capture_count",
                "fn": P.PARAMS["analysis_capture_oligos"],
                "action": "count",
                "overlap_fraction": 0.9,
            },
            {
                "name": "blacklist",
                "fn": P.PARAMS["analysis_blacklist"],
                "action": "count",
                "overlap_fraction": 1e-9,
            },
        ]
    ),
    r"ccanalysis/annotations/\1.annotations.tsv.gz",
)
def annotate_slices(infile, outfile):

    a = infile[0]
    colnames, fnames, actions, fractions = zip(
        *[(d["name"], d["fn"], d["action"], d["overlap_fraction"]) for d in infile[1]]
    )

    colnames = " ".join(colnames)
    fnames = " ".join([fn if not is_none(fn) else "-" for fn in fnames])
    actions = " ".join(actions)
    fractions = " ".join([str(frac) for frac in fractions])

    statement = f"""ccanalyser ccanalysis annotate_slices
                    -a {a}
                    -b {fnames}
                    --actions {actions}
                    --colnames {colnames}
                    -f {fractions}
                    -o {outfile}"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=6,
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(
    annotate_slices, mkdir("ccanalysis/reporters"), mkdir("ccanalysis/stats"),
)
@transform(
    align_reads,
    regex(r"aligned/(.*).bam"),
    add_inputs(r"ccanalysis/annotations/\1.annotations.tsv.gz"),
    r"ccanalysis/reporters/\1.fragments.tsv.gz",
)
def ccanalyser(infiles, outfile):
    """Processes bam files and annotations, filteres slices and outputs
    reporter slices for each capture site"""

    bam, annotations = infiles
    output_prefix = outfile.replace(".fragments.tsv.gz", "")
    stats_prefix = outfile.replace("/captures_and_reporters/", "/stats/").replace(
        ".fragments.tsv.gz", ".slice.stats"
    )

    statement = """ccanalyser ccanalysis ccanalysis
                    -i %(bam)s
                    -a %(annotations)s
                    --output_prefix %(output_prefix)s
                    --stats_out %(stats_prefix)s
                    --method %(analysis_method)s
                    > %(output_prefix)s.log 2>&1"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_advanced_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(ccanalyser, mkdir("ccanalysis/reporters_deduplicated"))
@collate(
    "ccanalysis/reporters/*.tsv.gz",
    regex(r"ccanalysis/reporters/(.*)_part\d+.(?:flashed|pe).fragments.tsv.gz"),
    r"ccanalysis/reporters_deduplicated/\1.hdf5",
)
def deduplicate_fragments(infiles, outfile):

    infiles = " ".join(infiles)
    statement = """ccanalyser ccanalysis rmdup_slices
                   fragments
                   -i %(infiles)s
                   -f %(outfile)s
                   --shuffle
                   -p %(pipeline_n_cores)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory="32G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(deduplicate_fragments)
@transform(
    "ccanalysis/reporters/*.tsv.gz",
    regex(r"ccanalysis/reporters/(.*)_(part\d+).(flashed|pe).(?!fragments)(.*).tsv.gz"),
    add_inputs(r"ccanalysis/reporters_deduplicated/\1.hdf5"),
    r"ccanalysis/reporters_deduplicated/\1.\2.\3.\4.tsv.gz",
)
def deduplicate_slices(infile, outfile):

    tsv, dedup_frag_ids = infile
    statement = """ccanalyser ccanalysis rmdup_slices
                   slices
                   -i %(tsv)s
                   -f %(dedup_frag_ids)s
                   -o %(outfile)s
                   -p %(pipeline_n_cores)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory=P.PARAMS["pipeline_advanced_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("ccanalysis/reporters_combined"))
@collate(
    deduplicate_slices,
    regex(r"ccanalysis/reporters_deduplicated/(.*)_part\d+.(.*).tsv.gz"),
    r"ccanalysis/reporters_combined/\1.\2.tsv.gz",
)
def collate_ccanalyser(infiles, outfile):
    """Combines multiple capture site tsv files"""

    inlist = " ".join(infiles)

    statement = """ccanalyser utils agg_tsv
                   concatenate
                   -i %(inlist)s
                   -o %(outfile)s
                   --header"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


# @active_if(P.PARAMS["analysis_method"] in ["tri", "tiled"])
@follows(mkdir("ccanalysis/interaction_counts"))
@transform(
    deduplicate_slices,
    regex(r"ccanalysis/reporters_deduplicated/(.*).tsv.gz"),
    r"ccanalysis/interaction_counts/\1.tsv.gz",
)
def count_rf_interactions(infile, outfile):

    statement = [
        "ccanalyser ccanalysis count_interactions",
        f"-f {infile}",
        f"-o {outfile}",
        "--remove_exclusions",
        f"> {outfile}.log",
    ]

    statement = " ".join(statement)

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=2,
        job_condaenv=P.PARAMS["conda_env"],
    )


@collate(
    count_rf_interactions,
    regex(r"ccanalysis/interaction_counts/(.*).part\d+.(?:flashed|pe).(.*).tsv.gz"),
    r"ccanalysis/interaction_counts_combined/\1.\2.tsv.gz",
)
def collate_rf_counts(infiles, outfile):
    """Combines multiple capture site tsv files"""

    inlist = " ".join(infiles)
    mkdir("ccanalysis/rf_counts_aggregated")

    statement = """ccanalyser utils agg_tsv
                   concatenate
                   -i %(inlist)s
                   --header
                   -o %(outfile)s.tmp.tsv
                   &&
                   ccanalyser utils agg_tsv
                   aggregate
                   -i %(outfile)s.tmp.tsv
                   --header
                   --groupby_columns rf1 rf2
                   --aggregate_columns count
                   --aggregate_method sum
                   -o %(outfile)s 
                   &&
                   rm %(outfile)s.tmp.tsv
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(
    collate_rf_counts,
    regex(r"ccanalysis/interaction_counts_combined/(.*).tsv.gz"),
    add_inputs(digest_genome),
    r"ccanalysis/interaction_counts_combined/\1_restriction_fragments.hdf5",
)
def store_rf_counts(infile, outfile):

    rf_counts, rf_map = infile
    outfile_prefix = outfile.replace("_restriction_fragments.hdf5", "")
    bin_sizes = " ".join(re.split(r"[,;]\s*|\s+", str(P.PARAMS["plot_bin_size"])))
    bin_method = P.PARAMS.get("plot_bin_method", "overlap")

    statement = """ccanalyser ccanalysis store_interactions
                   -c %(rf_counts)s
                   -m %(rf_map)s
                   -g %(genome_chrom_sizes)s
                   -o %(outfile_prefix)s
                   -b %(bin_sizes)s
                   --bin_method %(bin_method)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(store_rf_counts, mkdir("visualise"), mkdir("visualise/interaction_plots"))
@transform(
    "ccanalysis/interaction_counts_combined/*.hdf5",
    regex(r"ccanalysis/interaction_counts_combined/(.*)_(.*)_binned.hdf5"),
    add_inputs(r"ccanalysis/interaction_counts_combined/\1_restriction_fragments.hdf5"),
    r"visualise/interaction_plots/\1.log",
)
def plot_rf_counts(infile, outfile):

    if not is_off(P.PARAMS["plot_coordinates"]):

        cooler_binned, cooler_rf = infile
        output_prefix = outfile.replace(".log", "")
        normalisation = (
            P.PARAMS["plot_advanced_normalisation"]
            if not is_none(P.PARAMS["plot_advanced_normalisation"])
            else "infer"
        )
        normalisation_options = (
            P.PARAMS["plot_advanced_normalisation_options"]
            if not is_none(P.PARAMS["plot_advanced_normalisation_options"])
            else ""
        )

        statement = f"""ccanalyser ccanalysis plot_interactions
                       {cooler_rf}
                       {cooler_binned}
                       -c %(plot_coordinates)s
                       --method %(analysis_method)s
                       --normalisation %(normalisation)s
                       %(normalisation_options)s
                       --cmap %(plot_cmap)s
                       -o %(output_prefix)s
                       -f %(plot_format)s
                       > %(outfile)s 2>&1"""

        P.run(
            statement,
            job_queue=P.PARAMS["pipeline_cluster_queue"],
            job_condaenv=P.PARAMS["conda_env"],
        )

    else:
        print("Not plotting as no coordinates provided")


@active_if(P.PARAMS["analysis_method"] in ["capture", "tri"])
@follows(mkdir("ccanalysis/bedgraphs"))
@transform(
    collate_ccanalyser,
    regex(r"ccanalysis/reporters_combined/(.*).tsv.gz"),
    add_inputs(digest_genome),
    r"ccanalysis/bedgraphs/\1.bedgraph.gz",
)
def make_bedgraph(infiles, outfile):
    """Intersect reporters with genome restriction fragments to create bedgraph"""
    tsv_fn = infiles[0]
    re_map = infiles[1]
    statement = """ccanalyser ccanalysis slices_to_bdg
                   -i %(tsv_fn)s
                   -b %(re_map)s
                   -o %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("visualise/bigwigs"))
@transform(
    make_bedgraph,
    regex(r"ccanalysis/bedgraphs/(.*).bedgraph.gz"),
    r"visualise/bigwigs/\1.bigWig",
)
def make_bigwig(infile, outfile):
    """Uses UCSC tools bedGraphToBigWig to generate bigWigs for each bedgraph"""

    tmp = infile.replace(".gz", "")

    statement = """zcat %(infile)s > %(tmp)s
                   && bedGraphToBigWig %(tmp)s %(genome_chrom_sizes)s %(outfile)s
                   && rm %(tmp)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    if os.path.exists("chrom_sizes.txt.tmp"):
        os.remove("chrom_sizes.txt.tmp")


# @follows(ccanalyser)
# @merge(
#     [
#         "fastq/deduplicated/*.stats",
#         "fastq/digested/*.stats",
#         "ccanalysis/stats/*.slice.stats",
#         "ccanalysis/stats/*.reporter.stats",
#     ],
#     "run_statistics/combined_stats.tsv",
# )
# def aggregate_stats(infiles, outfile):

#     dedup = " ".join([fn for fn in infiles if "deduplicated/" in fn])
#     digestion = " ".join([fn for fn in infiles if "digested/" in fn])
#     slices = " ".join([fn for fn in infiles if ".slice.stats" in fn])
#     reporters = " ".join([fn for fn in infiles if ".reporter.stats" in fn])
#     outdir = os.path.dirname(outfile)

#     statement = """ccanalyser stats agg_stats
#                     --deduplication_stats %(dedup)s
#                     --digestion_stats %(digestion)s
#                     --ccanalyser_stats %(slices)s
#                     --reporter_stats %(reporters)s
#                     --output_dir %(outdir)s
#                 """

#     P.run(
#         statement,
#         job_queue=P.PARAMS["pipeline_cluster_queue"],
#         job_condaenv=P.PARAMS["conda_env"],
#     )


# @follows(aggregate_stats)
# @merge("run_statistics/*.tsv", r"run_statistics/visualise_run_statistics.html")
# def build_report(infile, outfile):
#     """Run jupyter notebook for reporting and plotting. First moves the notebook
#     then converts to html"""

#     statement = """rm run_statistics/visualise_run_statistics* -f &&
#                    papermill
#                    %(PACKAGE_DIR)s/stats/visualise_capture-c_stats.ipynb
#                    run_statistics/visualise_run_statistics.ipynb
#                    -p directory $(pwd)/run_statistics/ &&
#                    jupyter nbconvert
#                    --no-input
#                    --to html
#                    run_statistics/visualise_run_statistics.ipynb
#                    run_statistics/visualise_run_statistics.html
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
#                 "visualise_run_statistics.html",
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


def main(argv=None):
    if argv is None:
        argv = sys.argv

    P.main(argv)


if __name__ == "__main__":
    main()
