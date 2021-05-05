#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This pipeline processes data from Capture-C/NG Capture-C/Tri-C and Tiled-C sequencing
protocols designed to identify 3D interactions in the genome from a specified viewpoint.

It takes Illumina paired-end sequencing reads in fastq format 
(gzip compression is prefered) as input and performs the following steps:

1. Identifies all restriction fragments in the genome
2. Quality control of raw reads (fastqc, multiqc)
3. Splits fastqs into smaller files to enable fast parallel processing.
4. Removal of PCR duplicates based on exact sequence matches from fastq files
5. Trimming of reads to remove adaptor sequence (trim_galore)
6. Combining overlapping read pairs (FLASh)
7. In silico digestion of reads in fastq files
8. Alignment of fastq files with a user specified aligner (i.e. bowtie/bowtie2; BWA is not supported)
9. Analysis of alignment statistics (picard CollectAlignmentSummaryMetrics, multiqc)
10. Annotation of mapped reads with overlaps of capture probes, exclusion regions, blacklist, restriction fragments
11. Removal of non-reporter slices and indentification of reporters
12. Removal of PCR duplicates (exact coordinate matches) 
13. Storage of reporters in `cooler format <https.//cooler.readthedocs.io/en/latest/datamodel.html>`
14. Generation of bedgraphs/BigWigs.
15. Collation of run statistics and generation of a run report


Optional:

* Generation of a UCSC track hub for visualisation.
* Differential interaction identification.
* Generation of subtraction bedgraphs for between condition comparisons
* Plotting of heatmaps.  


@authors: asmith, dsims
"""

from collections import defaultdict
import os
import re
import sys
import pickle
from cgatcore import pipeline as P
from cgatcore.iotools import touch_file, zap_file
import itertools
import warnings
import glob
from cgatcore.pipeline.parameters import PARAMS

warnings.simplefilter("ignore", category=RuntimeWarning)

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
    originate,
    split,
)
from ccanalyser.tools.statistics import (
    collate_slice_data,
    collate_read_data,
    collate_cis_trans_data,
    collate_histogram_data,
    extract_trimming_stats,
)

from ccanalyser.utils import is_on, is_none, is_valid_bed, make_group_track


##############################
#   Set-up global parameters #
##############################

# Set up global parameters dict
P.get_parameters("config.yml")

# Determine the number of unique sample names
N_SAMPLES = len(
    {re.match(r"(.*)_R*[12].fastq.*", fn).group(1) for fn in glob.glob("*.fastq*")}
)

# Has valid plotting coordinate bed file
VALID_PLOT_COORDINATES = is_valid_bed(P.PARAMS.get("plot_coordinates"), verbose=False)

# Create a UCSC hub or not
MAKE_HUB = is_on(P.PARAMS.get("hub_create"))


##############################
#  Pipeline set-up functions #
##############################


def modify_pipeline_params_dict():

    """
    Modifies P.PARAMS dictionary.

    * Selects the correct conda enviroment
    * Ensures the correct queue manager is selected.
    * Corrects the name of a UCSC hub by removing spaces and incorrect characters.

    """

    # Fix cgat-core bugs
    P.PARAMS["cluster_queue_manager"] = P.PARAMS.get("pipeline_cluster_queue_manager")
    P.PARAMS["conda_env"] = P.PARAMS.get(
        "conda_env", os.path.basename(os.environ["CONDA_PREFIX"])
    )

    # Sanitise hub name
    if P.PARAMS["hub_name"]:
        P.PARAMS["hub_name"] = re.sub(r"[,\s+\t;:]", "_", P.PARAMS["hub_name"])

    # Convert entries to the correct python type
    for key in P.PARAMS:
        if is_none(P.PARAMS[key]):
            P.PARAMS[key] = None
        elif is_on(P.PARAMS):
            P.PARAMS[key] = True


def set_up_chromsizes():
    """
    Ensures that genome chromsizes are present.

    If chromsizes are not provided this function attempts to download them from UCSC.
    The P.PARAMS dictionary is updated with the location of the chromsizes.

    """

    assert P.PARAMS.get("genome_name"), "Genome name has not been provided."

    if P.PARAMS["genome_chrom_sizes"]:
        pass

    elif os.path.exists("chrom_sizes.txt.tmp"):
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"

    else:
        from pybedtools.helpers import get_chromsizes_from_ucsc

        get_chromsizes_from_ucsc(P.PARAMS["genome_name"], "chrom_sizes.txt.tmp")
        P.PARAMS["genome_chrom_sizes"] = "chrom_sizes.txt.tmp"


def check_user_supplied_paths():

    paths_to_check = [
        "genome_fasta",
        "genome_aligner_index",
        "analysis_viewpoints",
    ]

    chrom_sizes = P.PARAMS["genome_chrom_sizes"]
    if any(ext in chrom_sizes for ext in [".txt", ".fai", ".tsv"]):
        paths_to_check.append("genome_chrom_sizes")

    if MAKE_HUB:
        paths_to_check.append("hub_dir")

    for path_name in paths_to_check:

        path_supplied = P.PARAMS[path_name]

        if not os.path.exists(path_supplied):

            if path_name == "genome_aligner_index":
                indicies = glob.glob(path_supplied + "*")
                if not len(indicies) >= 1:
                    raise OSError(f"Supplied indicies at: {path_supplied} do not exist")

            else:
                raise OSError(
                    f"Supplied path for {path_name}: {path_supplied} does not exist"
                )


##################
# Prepare genome #
#################


@follows(mkdir("ccanalyser_preprocessing/restriction_enzyme_map/"))
@transform(
    P.PARAMS.get("genome_fasta"),
    regex(r".*/(.*).fa.*"),
    r"ccanalyser_preprocessing/restriction_enzyme_map/genome.digest.bed.gz",
)
def genome_digest(infile, outfile):
    """
    In silco digestion of the genome to identify restriction fragment coordinates.

    Runs :ref:`ccanalyser genome digest <CLI Documentation>`.

    """

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


#############################
# Fastq file pre-processing #
#############################


@follows(mkdir("ccanalyser_preprocessing"), mkdir("ccanalyser_preprocessing/fastqc"))
@transform(
    "*.fastq*", regex(r"(.*).fastq.*"), r"ccanalyser_preprocessing/fastqc/\1_fastqc.zip"
)
def fastq_qc(infile, outfile):
    """Runs fastqc on the input files to generate fastq statistics."""

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
@merge(fastq_qc, "statistics/fastqc_report.html")
def fastq_multiqc(infile, outfile):
    """Collate fastqc reports into single report using multiqc"""

    bn = os.path.basename(outfile)
    dn = os.path.dirname(outfile)

    s1 = "rm -f %(outfile)s &&"
    s2 = "export LC_ALL=en_US.UTF-8 &&"
    s3 = "export LANG=en_US.UTF-8 &&"
    s4 = " ".join(
        ["multiqc", "ccanalyser_preprocessing/fastqc/", "-o", "%(dn)s", "-n", "%(bn)s"]
    )

    statement = " ".join([s1, s2, s3, s4])

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("ccanalyser_preprocessing/split"))
@collate(
    "*.fastq.gz",
    regex(r"(.*)_R*[12].fastq.*"),
    r"ccanalyser_preprocessing/split/\1.completed",
)
def fastq_split(infiles, outfile):
    """
    Splits the input fastq files into chunks for parallel processing

    Runs :ref:`ccanalyser fastq split <CLI Documentation>`.

    """

    infiles = " ".join(infiles)
    output_prefix = outfile.replace(".log", "")

    statement = """ccanalyser fastq split
                %(infiles)s
                -m unix
                -o %(output_prefix)s
                -n %(split_n_reads)s
                --no-gzip
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Create sentinel file
    touch_file(outfile)


@active_if(P.PARAMS.get("deduplication_pre-dedup", False))
@follows(
    mkdir("ccanalyser_preprocessing/deduplicated"),
    mkdir("ccanalyser_preprocessing/deduplicated/deduplicated_ids"),
    fastq_split,
)
@collate(
    "ccanalyser_preprocessing/split/*.fastq*",
    regex(r"ccanalyser_preprocessing/split/(.*)_part(\d+)_[12].fastq(?:.gz)?"),
    r"ccanalyser_preprocessing/deduplicated/deduplicated_ids/\1_\2.json.gz",
    extras=[r"\1", r"\2"],
)
def fastq_duplicates_parse(infiles, outfile, sample_name, part_no):

    """
    Parses fastq files into json format for sequence deduplication.

    Runs :ref:`ccanalyser fastq deduplicate parse <CLI Documentation>`

    """

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
    regex(r"ccanalyser_preprocessing/deduplicated/deduplicated_ids/(.*)_\d*.json.gz"),
    r"ccanalyser_preprocessing/deduplicated/deduplicated_ids/\1.json.gz",
)
def fastq_duplicates_identify(infiles, outfile):

    """
    Identifies duplicate sequences from parsed fastq files in json format.

    Runs :ref:`ccanalyser fastq deduplicate identify <CLI Documentation>`

    """

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


@active_if(P.PARAMS.get("deduplication_pre-dedup", False))
@follows(
    fastq_duplicates_parse,
    fastq_duplicates_identify,
    mkdir("statistics/deduplication/data/"),
)
@collate(
    "ccanalyser_preprocessing/split/*.fastq*",
    regex(r".*/(.*_part\d+)_[12].fastq(?:.gz)?"),
    r"ccanalyser_preprocessing/deduplicated/\1.completed",
)
def fastq_duplicates_remove(infiles, outfile):

    """
    Removes duplicate read fragments identified from parsed fastq files.
    """

    fq1, fq2 = infiles
    sample = re.match(r".*/(.*)(_part\d+)_[12].fastq(?:.gz)?", fq1)
    sample_name = sample.group(1)
    sample_part = sample.group(2)
    dd_ids = (
        f"ccanalyser_preprocessing/deduplicated/deduplicated_ids/{sample_name}.json.gz"
    )
    stats_prefix = f"statistics/deduplication/data/{sample_name}_{sample_part}"
    output_prefix = outfile.replace(".completed", "")

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

    # Make sentinel file
    touch_file(outfile)

    # Replace infiles with empty files
    for fn in infiles:
        zap_file(fn)


@follows(fastq_duplicates_remove)
@merge(
    "statistics/deduplication/data/*.csv",
    "statistics/deduplication/deduplication.reads.csv",
)
def stats_deduplication_collate(infiles, outfile):

    """Combines deduplication statistics from fastq file partitions."""

    stats_prefix = outfile.replace(".reads.csv", "")

    df_stats = collate_read_data(infiles)

    df_stats_read = df_stats.query('stat_type != "reads_removed"')

    df_stats.to_csv(f"{stats_prefix}.summary.csv", index=False)
    df_stats_read.to_csv(
        outfile, index=False
    )  # Modified to enable more streamlined summary at final stage


@follows(
    mkdir("ccanalyser_preprocessing/trimmed"),
    fastq_duplicates_remove,
    mkdir("statistics/trimming/data/"),
)
@collate(
    "ccanalyser_preprocessing/deduplicated/*.fastq*",
    regex(r"ccanalyser_preprocessing/deduplicated/(.*)_[12].fastq(?:.gz)?"),
    r"ccanalyser_preprocessing/trimmed/\1.completed",
)
def fastq_trim(infiles, outfile):

    """Trim adaptor sequences from fastq files using trim_galore"""

    fq1, fq2 = infiles
    fq1_basename, fq2_basename = os.path.basename(fq1), os.path.basename(fq2)

    outdir = os.path.dirname(outfile)
    trim_options = P.PARAMS.get("trim_options", " ")

    statement = """trim_galore
                   --cores %(pipeline_n_cores)s
                   --paired %(trim_options)s
                   --gzip
                   -o %(outdir)s
                   %(fq1)s
                   %(fq2)s
                   && mv ccanalyser_preprocessing/trimmed/%(fq1_basename)s_trimming_report.txt statistics/trimming/data
                   && mv ccanalyser_preprocessing/trimmed/%(fq2_basename)s_trimming_report.txt statistics/trimming/data
                   """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Make sentinel file
    touch_file(outfile)

    # Zero input files to save space
    for fn in infiles:
        zap_file(fn)


@follows(fastq_trim)
@merge("statistics/trimming/data/*.txt", r"statistics/trimming/trimming.summary.csv")
def stats_trim_collate(infiles, outfile):

    """Extracts and collates adapter trimming statistics from trim_galore output"""

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


@follows(fastq_trim, mkdir("ccanalyser_preprocessing/flashed"))
@collate(
    "ccanalyser_preprocessing/trimmed/*.fq*",
    regex(r"ccanalyser_preprocessing/trimmed/(.*)_[12]_.*.fq(?:.gz)?"),
    r"ccanalyser_preprocessing/flashed/\1.completed",
)
def fastq_flash(infiles, outfile):

    """Combine overlapping paired-end reads using FLASh"""

    fq1, fq2 = infiles
    output_prefix = outfile.replace(".completed", "")
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

    touch_file(outfile)


@follows(
    mkdir("ccanalyser_preprocessing/digested"),
    fastq_flash,
    mkdir("statistics/digestion/data"),
)
@transform(
    "ccanalyser_preprocessing/flashed/*.fastq.gz",
    regex(r"ccanalyser_preprocessing/flashed/(.*).extendedFrags.fastq.gz"),
    r"ccanalyser_preprocessing/digested/\1.flashed.fastq.gz",
)
def fastq_digest_combined(infile, outfile):

    """In silico restriction enzyme digest of combined (flashed) read pairs"""

    fn = os.path.basename(infile).replace(".extendedFrags.fastq.gz", "")
    sn = fn.split("_part")[0]
    n_cores = P.PARAMS["pipeline_n_cores"]
    n_cores_digestion = int(n_cores) - 2 if int(n_cores) > 2 else 1

    statement = """ccanalyser fastq digest
                   %(infile)s
                   -m flashed
                   -r %(analysis_restriction_enzyme)s
                   -o %(outfile)s
                   --stats_prefix %(outfile)s.stats
                   --minimum_slice_length 18
                   --compression_level %(pipeline_compression)s
                   -p 1
                   --stats_prefix statistics/digestion/data/%(fn)s.flashed
                   --sample_name %(sn)s
                    """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=3,
        job_condaenv=P.PARAMS["conda_env"],
    )

    zap_file(infile)


@follows(fastq_flash, mkdir("statistics/digestion"))
@collate(
    "ccanalyser_preprocessing/flashed/*.fastq.gz",
    regex(r"ccanalyser_preprocessing/flashed/(.*).notCombined_[12].fastq.gz"),
    r"ccanalyser_preprocessing/digested/\1.pe.fastq.gz",
)
def fastq_digest_non_combined(infiles, outfile):

    """In silico restriction enzyme digest of non-combined (non-flashed) read pairs"""

    fq1, fq2 = infiles
    fn = os.path.basename(fq1).replace(".notCombined_1.fastq.gz", "")
    sn = fn.split("_part")[0]
    n_cores = P.PARAMS["pipeline_n_cores"]
    n_cores_digestion = int(n_cores) - 2 if int(n_cores) > 2 else 1

    statement = """ccanalyser fastq digest
                   %(fq1)s
                   %(fq2)s
                   -m pe
                   -r %(analysis_restriction_enzyme)s
                   -o %(outfile)s
                   --minimum_slice_length 18
                   -p 1
                   --compression_level %(pipeline_compression)s
                   --stats_prefix statistics/digestion/data/%(fn)s.pe
                   --sample_name %(sn)s
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=3,
        job_condaenv=P.PARAMS["conda_env"],
    )

    for fn in infiles:
        zap_file(fn)


@follows(fastq_digest_combined, fastq_digest_non_combined)
@merge("statistics/digestion/data/*", "statistics/digestion/digestion.reads.csv")
def stats_digestion_collate(infiles, outfile):

    """Aggregates in silico digestion statistics from fastq file partitions."""

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

    # Collate histogram, read and slice statistics
    df_hist_filt = collate_histogram_data(data["hist_filt"])
    df_hist_unfilt = collate_histogram_data(data["hist_unfilt"])
    # df_slice = collate_read_data(data["slice"])
    df_read = collate_read_data(data["read"])

    # Merge filtered and unfiltered histograms
    df_hist = pd.concat(
        [df_hist_unfilt.assign(filtered=0), df_hist_filt.assign(filtered=1)]
    ).sort_values(["sample", "read_type", "n_slices"])

    # Output histogram, slice and read statics
    df_hist.to_csv(f"{stats_prefix}.histogram.csv", index=False)
    # df_slice.to_csv(f"{stats_prefix}.slice.csv", index=False)
    df_read.to_csv(outfile, index=False)


@follows(fastq_digest_combined, fastq_digest_non_combined)
def fastq_preprocessing():
    pass


#################################
# Read alignment and processing #
#################################


@follows(mkdir("ccanalyser_preprocessing"), fastq_preprocessing)
@transform(
    [fastq_digest_combined, fastq_digest_non_combined],
    regex(r"ccanalyser_preprocessing/digested/(.*).fastq.gz"),
    r"ccanalyser_preprocessing/aligned/\1.bam",
)
def fastq_alignment(infile, outfile):

    """Aligns in silico digested fastq files to the genome."""

    aligner = P.PARAMS["align_aligner"]
    index_flag = P.PARAMS.get("align_index_flag", " ")
    options = P.PARAMS.get("align_options", " ")

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
    fastq_alignment,
    regex(r"ccanalyser_preprocessing/aligned/(.*)_part\d+.*.bam"),
    r"ccanalyser_preprocessing/aligned/\1.bam",
)
def alignments_merge(infiles, outfile):
    """
    Combines bam files (by flashed/non-flashed status and sample).

    This task simply provides an input for picard CollectAlignmentSummaryMetrics
    and is only used to provide overall mapping statistics. Fastq partitions
    are *not* combined at this stage.

    """

    fnames = " ".join(infiles)

    statement = """samtools merge %(outfile)s %(fnames)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(alignments_merge, regex("(.*).bam"), r"\1.bam.bai")
def alignments_index(infile, outfile):

    """Indexes all bam files (both partitioned and merged)"""

    statement = """samtools index %(infile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="1G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(fastq_alignment)
def pre_annotation():
    pass


############################
# Annotation of alignments #
############################


@originate("ccanalyser_analysis/annotations/exclude.bed")
def annotate_make_exclusion_bed(outfile):

    """Generates exclusion window around each capture site"""

    assert is_valid_bed(P.PARAMS['analysis_viewpoints']), "Viewpoints bed file is not a valid bed file"
    statement = """bedtools slop
                    -i %(analysis_viewpoints)s -g %(genome_chrom_sizes)s -b %(analysis_reporter_exclusion_zone)s
                    | bedtools subtract -a - -b %(analysis_viewpoints)s
                    | sort -k1,1 -k2,2n 
                    > %(outfile)s"""
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@originate("ccanalyser_analysis/annotations/viewpoints.bed")
def annotate_sort_viewpoints(outfile):

    """Sorts the capture oligos for bedtools intersect with --sorted option"""
    
    assert is_valid_bed(P.PARAMS['analysis_viewpoints']), "Viewpoints bed file is not a valid bed file"
    statement = """cat %(analysis_viewpoints)s | sort -k1,1 -k2,2n > %(outfile)s"""
    
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(P.PARAMS.get("analysis_optional_blacklist"))
@originate("ccanalyser_analysis/annotations/blacklist.bed")
def annotate_sort_blacklist(outfile):

    """Sorts the capture oligos for bedtools intersect with --sorted option"""

    if os.path.exists(P.PARAMS["analysis_optional_blacklist"]):
        statement = (
            """cat %(analysis_optional_blacklist)s | sort -k1,1 -k2,2n > %(outfile)s"""
        )
    else:
        # Make a blank file if no blacklist
        statement = "touch %(outfile)s"

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(
    genome_digest,
    annotate_make_exclusion_bed,
    annotate_sort_viewpoints,
    annotate_sort_blacklist,
)
@transform(
    fastq_alignment,
    regex(r".*/(.*).bam"),
    add_inputs(
        [
            {
                "name": "restriction_fragment",
                "fn": "ccanalyser_preprocessing/restriction_enzyme_map/genome.digest.bed.gz",
                "action": "get",
                "fraction": 0.2,
            },
            {
                "name": "capture",
                "fn": "ccanalyser_analysis/annotations/viewpoints.bed",
                "action": "get",
                "fraction": 0.9,
            },
            {
                "name": "exclusion",
                "fn": "ccanalyser_analysis/annotations/exclude.bed",
                "action": "get",
                "fraction": 1e-9,
            },
            {
                "name": "exclusion_count",
                "fn": "ccanalyser_analysis/annotations/exclude.bed",
                "action": "count",
                "fraction": 1e-9,
            },
            {
                "name": "capture_count",
                "fn": "ccanalyser_analysis/annotations/viewpoints.bed",
                "action": "count",
                "fraction": 0.9,
            },
            {
                "name": "blacklist",
                "fn": "ccanalyser_analysis/annotations/blacklist.bed",
                "action": "count",
                "fraction": 1e-9,
            },
        ]
    ),
    r"ccanalyser_analysis/annotations/\1.annotations.tsv",
)
def annotate_alignments(infile, outfile):

    """
    Annotates mapped read slices.

    Slices are annotated with:
     * capture name
     * capture count
     * exclusion name
     * exclusion count
     * blacklist count
     * restriction fragment number
    """

    slices = infile[0]
    flags = {"name": "-n", "fn": "-b", "action": "-a", "fraction": "-f"}

    cmd_args = []
    for args in infile[1]: # infile[1] == arguments for annotation
        for arg_name, arg in args.items():
            cmd_args.append(f'{flags.get(arg_name)} {arg if arg else "-"}')

    cmd_args = " ".join(cmd_args)
    statement = """bedtools bamtobed -i %(slices)s
                   | sort -k1,1 -k2,2n 
                   | ccanalyser alignments annotate
                   -
                   %(cmd_args)s
                   -o %(outfile)s
                   --invalid_bed_action ignore
                   -p 1
                """

    P.run(
        statement.replace("\n", " "),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=1,
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(fastq_preprocessing, annotate_alignments)
def post_annotation():
    """Runs the pipeline until just prior to identification of reporters"""
    pass


###########################
# Filtering of alignments #
###########################


@follows(
    post_annotation,
    annotate_alignments,
    mkdir("ccanalyser_analysis/reporters/unfiltered/"),
    mkdir("statistics/reporters/data"),
)
@transform(
    fastq_alignment,
    regex(r"ccanalyser_preprocessing/aligned/(.*).bam"),
    add_inputs(r"ccanalyser_analysis/annotations/\1.annotations.tsv"),
    r"ccanalyser_analysis/reporters/unfiltered/\1.completed",
)
def alignments_filter(infiles, outfile):
    """Filteres slices and outputs reporter slices for each capture site"""

    bam, annotations = infiles
    sample = re.match(r".*/(.*)_(part\d+).(flashed|pe).bam", bam)
    sample_name = sample.group(1)
    sample_part = sample.group(2)
    sample_read_type = sample.group(3)

    output_prefix = outfile.replace(".completed", "")
    output_log_file = f'{output_prefix}.log'
    stats_prefix = (
        f"statistics/reporters/data/{sample_name}_{sample_part}_{sample_read_type}"
    )

    statement = """ccanalyser
                   alignments
                   filter
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
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Make sentinel file
    touch_file(outfile)

    # Zero annotations
    if not P.PARAMS.get('analysis_optional_keep_annotations', False):
        zap_file(annotations)


@follows(mkdir("ccanalyser_analysis/reporters/collated"), alignments_filter)
@collate(
    "ccanalyser_analysis/reporters/unfiltered/*.tsv",
    regex(
        r".*/(?P<sample>.*)_part\d+.(flashed|pe).(?P<capture>.*).(slices|fragments).tsv"
    ),
    r"ccanalyser_analysis/reporters/collated/\1.\2.\3.\4.tsv",
    extras=[r"\1", r"\2", r"\3", r"\4"],
)
def reporters_collate(infiles, outfile, *grouping_args):

    """Concatenates identified reporters """

    statement = []
    for ii, fn in enumerate(infiles):
        if ii == 0:
            cmd = f"cat {fn} > {outfile}"
        else:
            cmd = f"tail -n +2 {fn} >> {outfile}"

        statement.append(cmd)

    P.run(
        " && ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zero un-aggregated reporters
    for fn in infiles:
        zap_file(fn)


@follows(alignments_filter, mkdir("ccanalyser_analysis/reporters/deduplicated"))
@transform(
    reporters_collate,
    regex(r".*/(?P<sample>.*).(flashed|pe).(?P<capture>.*).fragments.tsv"),
    r"ccanalyser_analysis/reporters/deduplicated/\1.\2.\3.json.gz",
    extras=[r"\2"],
)
def alignments_deduplicate_fragments(infile, outfile, read_type):

    """
    Identifies duplicate fragments with the same coordinates and order.
    """

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


@follows(alignments_deduplicate_fragments)
@transform(
    reporters_collate,
    regex(
        r"ccanalyser_analysis/reporters/collated/(.*)\.(flashed|pe)\.(.*)\.slices.tsv"
    ),
    add_inputs(r"ccanalyser_analysis/reporters/deduplicated/\1.\2.\3.json.gz"),
    r"ccanalyser_analysis/reporters/deduplicated/\1.\2.\3.slices.tsv",
    extras=[r"\1", r"\2", r"\3"],
)
def alignments_deduplicate_slices(
    infile, outfile, sample_name, read_type, capture_oligo
):

    """Removes reporters with duplicate coordinates"""

    slices, duplicated_ids = infile
    stats_prefix = (
        f"statistics/reporters/data/{sample_name}_{read_type}_{capture_oligo}"
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
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Zero non-deduplicated reporters
    zap_file(slices)


@collate(
    alignments_deduplicate_slices,
    regex(r".*/(?P<sample>.*).(?:flashed|pe).(?P<capture>.*).slices.tsv"),
    r"ccanalyser_analysis/reporters/\1.\2.tsv.gz",
    extras=[r"\1", r"\2"],
)
def alignments_deduplicate_collate(infiles, outfile, *grouping_args):

    """Final collation of reporters by sample and capture probe"""

    statement = []
    tmp = outfile.replace(".gz", "")
    for ii, fn in enumerate(infiles):
        if ii == 0:
            cmd = f"cat {fn} > {tmp}"
        else:
            cmd = f"tail -n +2 {fn} >> {tmp}"

        statement.append(cmd)

    statement.append(f'cat {tmp} | pigz -p {P.PARAMS["pipeline_n_cores"]} > {outfile}')

    P.run(
        " && ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(alignments_deduplicate_collate)
@merge("statistics/reporters/data/*", "statistics/reporters/reporters.reads.csv")
def stats_alignment_filtering_collate(infiles, outfile):

    """'Combination of all reporter identification and filtering statistics"""

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


@follows(alignments_deduplicate_slices, stats_alignment_filtering_collate)
def post_ccanalyser_analysis():
    """Reporters have been identified, deduplicated and collated by sample/capture probe"""


####################
# Reporter storage #
####################


@follows(mkdir("ccanalyser_analysis/reporters/counts"))
@transform(
    alignments_deduplicate_collate,
    regex(r"ccanalyser_analysis/reporters/(.*)\.(.*).tsv.gz"),
    r"ccanalyser_analysis/reporters/counts/\1.\2.tsv.gz",
)
def reporters_count(infile, outfile):

    """Counts the number of interactions identified between reporter restriction fragments"""

    statement = [
        "ccanalyser reporters  count",
        "%(infile)s",
        "-o %(outfile)s",
        "--remove_exclusions",
        "> %(outfile)s.log",
    ]

    statement = " ".join(statement)

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=2,
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("ccanalyser_analysis/reporters/fragments"))
@transform(
    reporters_count,
    regex(r"ccanalyser_analysis/reporters/counts/(.*)\.(.*)\.tsv.gz"),
    add_inputs(genome_digest),
    r"ccanalyser_analysis/reporters/fragments/\1.\2.fragments.hdf5",
    extras=[r"\1", r"\2"],
)
def reporters_store_restriction_fragment(infile, outfile, sample_name, capture_name):

    """Stores restriction fragment interaction counts in cooler format"""

    counts, rf_map = infile
    output_prefix = outfile.replace(f".{capture_name}.fragments", "")

    statement = """ccanalyser reporters  store
                   fragments
                   %(counts)s
                   -f %(rf_map)s
                   -g %(genome_name)s
                   -n %(capture_name)s
                   -o %(output_prefix)s
                   -c %(analysis_viewpoints)s
                   --suffix fragments
                   """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(genome_digest, reporters_count)
@originate(r"ccanalyser_analysis/reporters/binners.pkl")
def generate_bin_conversion_tables(outfile):
    """
    Converts restriction fragments to genomic bins.

    Binning restriction fragments into genomic bins takes a substantial
    amount of time and memory. To avoid repeatedly performing the same action,
    bin conversion tables are calculated once for each required resolution and
    then stored as a pickle file.

    """

    from ccanalyser.tools.storage import GenomicBinner

    frags = pd.read_csv(
        "ccanalyser_preprocessing/restriction_enzyme_map/genome.digest.bed.gz",
        sep="\t",
        names=["chrom", "start", "end", "name"],
    )

    binner_dict = dict()
    for bs in re.split(r"[,;]\s*|\s+", str(P.PARAMS["analysis_bin_size"])):
        gb = GenomicBinner(
            chromsizes=P.PARAMS["genome_chrom_sizes"], fragments=frags, binsize=int(bs)
        )
        gb.bin_conversion_table  # Property is cached so need to call it to make sure it is present.
        binner_dict[int(bs)] = gb

    with open("ccanalyser_analysis/reporters/binners.pkl", "wb") as w:
        pickle.dump(binner_dict, w)


@active_if(P.PARAMS.get("analysis_bin_size"))
@follows(generate_bin_conversion_tables, mkdir("ccanalyser_analysis/reporters/binned/"))
@transform(
    reporters_store_restriction_fragment,
    regex(r"ccanalyser_analysis/reporters/fragments/(.*)\.(.*)\.fragments\.hdf5"),
    add_inputs(generate_bin_conversion_tables),
    r"ccanalyser_analysis/reporters/binned/\1.\2.completed",
    extras=[r"\2"],
)
def reporters_store_binned(infile, outfile, capture_name):

    """
    Converts a cooler file of restriction fragments to even genomic bins.
    """

    infile, conversion_tables = infile

    bin_options = " -b " + " -b ".join(
        re.split(r"[,;]\s*|\s+", str(P.PARAMS["analysis_bin_size"]))
    )
    output_prefix = outfile.replace(f".{capture_name}.completed", "")

    statement = f"""ccanalyser reporters  store
                  bins
                  %(infile)s
                  -o %(output_prefix)s
                  %(bin_options)s
                  --conversion_tables %(conversion_tables)s
                  --normalise
                  -p 4
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=4,
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Make sentinel file
    touch_file(outfile)


@follows(reporters_store_restriction_fragment, reporters_store_binned)
@collate(
    [
        "ccanalyser_analysis/reporters/fragments/*.hdf5",
        "ccanalyser_analysis/reporters/binned/*.hdf5",
    ],
    regex(r".*/(.*)\.(.*)\.(?:fragments|\d+)\.hdf5"),
    r"ccanalyser_analysis/reporters/\1.hdf5",
    extras=[r"\1"],
)
def reporters_store_merged(infiles, outfile, sample_name):

    """Combines cooler files together"""

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


#######################
# Pipeline statistics #
#######################


@merge(
    [
        stats_deduplication_collate,
        stats_digestion_collate,
        stats_alignment_filtering_collate,
    ],
    "statistics/run_statistics.csv",
)
def pipeline_merge_stats(infiles, outfile):

    """
    Generates a summary statistics file for the pipeline run.

    """

    df = pd.concat([pd.read_csv(fn) for fn in infiles])

    df.sort_values(
        ["sample", "read_type", "stat"], ascending=[True, True, False]
    ).to_csv(outfile)


@merge(
    [pipeline_merge_stats],
    "statistics/visualise_statistics.html",
)
def pipeline_make_report(infile, outfile):
    """Run jupyter notebook for reporting and plotting pipeline statistics"""

    # Make sure black cache is generated
    import black

    black.CACHE_DIR.mkdir(parents=True, exist_ok=True)

    path_script = __file__
    path_script_dir = os.path.dirname(path_script)
    path_nb_dir = os.path.dirname(path_script_dir)

    statement_clean = "rm statistics/visualise_statistics* -f"

    statement_papermill = """papermill
                             -k python3
                             -p directory $(pwd)/statistics/
                             %(path_nb_dir)s/visualise_statistics.ipynb
                             statistics/visualise_statistics.ipynb
                            """
    statement_nbconvert = """jupyter nbconvert
                             --no-input
                             --to html
                             statistics/visualise_statistics.ipynb
                             statistics/visualise_statistics.html
                          """

    P.run(
        " && ".join([statement_clean, statement_papermill, statement_nbconvert]),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#####################
# Reporter pileups  #
#####################


@follows(mkdir("ccanalyser_analysis/bedgraphs"))
@transform(
    reporters_store_merged,
    regex(r".*/(.*).hdf5"),
    r"ccanalyser_analysis/bedgraphs/\1.raw.completed",
    extras=[r"\1"],
)
def reporters_make_bedgraph(infile, outfile, sample_name):
    """Extract reporters in bedgraph format from stored interactions"""

    output_prefix = f"ccanalyser_analysis/bedgraphs/{sample_name}.raw"

    statement = """ccanalyser reporters  pileup
                   %(infile)s
                   -o %(output_prefix)s
                """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    touch_file(outfile)


@transform(
    reporters_store_merged,
    regex(r".*/(.*).hdf5"),
    r"ccanalyser_analysis/bedgraphs/\1.normalised.completed",
    extras=[r"\1"],
)
def reporters_make_bedgraph_normalised(infile, outfile, sample_name):
    """
    Extract reporters in bedgraph format from stored interactions.

    In addition to generating a bedgraph this task also normalises the counts
    by the number of cis interactions identified to enable cross sample comparisons.

    """

    output_prefix = f"ccanalyser_analysis/bedgraphs/{sample_name}.normalised"

    statement = """ccanalyser reporters  pileup
                   %(infile)s
                   -o %(output_prefix)s
                   --normalise
                   """
    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    touch_file(outfile)


@active_if(N_SAMPLES >= 2)
@follows(
    mkdir("ccanalyser_compare/bedgraphs_union"),
    reporters_make_bedgraph,
    reporters_make_bedgraph_normalised,
)
@collate(
    "ccanalyser_analysis/bedgraphs/*.bedgraph",
    regex(r".*/(?:.*)\.(raw|normalised|windowed)\.(.*).bedgraph"),
    r"ccanalyser_compare/bedgraphs_union/\2.\1.tsv",
    extras=[r"\1", r"\2"],
)
def reporters_make_union_bedgraph(infiles, outfile, normalisation_type, capture_name):

    """
    Collates bedgraphs by capture probe into a single file for comparison.

    See `bedtools unionbedg <https://bedtools.readthedocs.io/en/latest/content/tools/unionbedg.html>`_
    for more details.

    """

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


@active_if(N_SAMPLES >= 2)
@follows(
    mkdir("ccanalyser_compare/bedgraphs_comparison/"), reporters_make_union_bedgraph
)
@transform(
    reporters_make_union_bedgraph,
    regex(r"ccanalyser_compare/bedgraphs_union/(.*)\.normalised\.tsv"),
    r"ccanalyser_compare/bedgraphs_comparison/\1.completed",
    extras=[r"\1"],
)
def reporters_make_comparison_bedgraph(infile, outfile, viewpoint):

    import numpy as np

    df_bdg = pd.read_csv(infile, sep="\t")
    dir_output = os.path.dirname(outfile)

    summary_methods = re.split(r'[,;\s+]', P.PARAMS.get('compare_summary_methods', 'mean'))
    summary_functions = {method: getattr(np, method) for method in summary_methods}

    # If no design matrix, make one assuming the format has been followed
    if not P.PARAMS.get("analysis_design"):
        col_dict = {col: "_".join(col.split("_")[:-1]) for col in df_bdg.columns[3:]}
        df_design = pd.Series(col_dict).to_frame("condition")

    condition_groups = df_design.groupby("condition").groups

    for a, b in itertools.permutations(condition_groups, 2):

        # Extract the two groups
        df_a = df_bdg.loc[:, condition_groups[a]]
        df_b = df_bdg.loc[:, condition_groups[b]]

        for summary_method in summary_functions:
            # Get summary counts
            a_summary = df_a.pipe(summary_functions[summary_method], axis=1)
            b_summary = df_b.pipe(summary_functions[summary_method], axis=1)

            df_a_bdg = pd.concat([df_bdg.iloc[:, :3], a_summary], axis=1)
            df_b_bdg = pd.concat([df_bdg.iloc[:, :3], b_summary], axis=1)
            df_subtraction_bdg = pd.concat([df_bdg.iloc[:, :3], a_summary - b_summary], axis=1)

            df_a_bdg.to_csv(f"{dir_output}/{a}_{summary_method}.{viewpoint}.bedgraph",
                            sep='\t',
                            header=False,
                            index=None)

            df_b_bdg.to_csv(f"{dir_output}/{b}_{summary_method}.{viewpoint}.bedgraph",
                            sep='\t',
                            header=False,
                            index=None)
            
            df_subtraction_bdg.to_csv(f"{dir_output}/{a}_vs_{b}.{summary_method}-subtraction.{viewpoint}.bedgraph",
                                      sep='\t',
                                      index=None,
                                      header=False)

    touch_file(outfile)


@follows(
    mkdir("ccanalyser_analysis/bigwigs"),
    reporters_make_bedgraph,
    reporters_make_bedgraph_normalised,
    reporters_make_comparison_bedgraph,
)
@transform(
    [
        "ccanalyser_analysis/bedgraphs/*",
        "ccanalyser_compare/bedgraphs_comparison/*.bedgraph",
    ],
    regex(r".*/(.*).bedgraph"),
    r"ccanalyser_analysis/bigwigs/\1.bigWig",
)
def reporters_make_bigwig(infile, outfile):
    """Uses UCSC tools bedGraphToBigWig to generate bigWigs for each bedgraph"""

    tmp = f"{outfile}.tmp"
    statement = """  cat %(infile)s
                   | sort -k1,1 -k2,2n > %(tmp)s
                   && bedGraphToBigWig %(tmp)s %(genome_chrom_sizes)s %(outfile)s
                   && rm %(tmp)s
                """

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#######################
# UCSC hub generation #
#######################


@active_if(MAKE_HUB)
@merge(
    reporters_make_bigwig,
    os.path.join(
        P.PARAMS.get("hub_dir", ""), P.PARAMS.get("hub_name", "") + ".hub.txt"
    ),
    extras=[pipeline_make_report],
)
def hub_make(infiles, outfile, statistics):
    """Creates a ucsc hub from the pipeline output"""

    import trackhub

    excluded = [
        "raw",
    ]

    bigwigs = [fn for fn in infiles if not any(e in fn for e in excluded)]
    key_sample = lambda b: os.path.basename(b).split(".")[0]
    key_capture = lambda b: b.split(".")[-2]

    # Need to make an assembly hub if this is a custom genome
    if not P.PARAMS.get("genome_custom"):

        hub, genomes_file, genome, trackdb = trackhub.default_hub(
            hub_name=P.PARAMS["hub_name"],
            short_label=P.PARAMS.get("hub_short"),
            long_label=P.PARAMS.get("hub_long"),
            email=P.PARAMS["hub_email"],
            genome=P.PARAMS["genome_name"],
        )

        for key in [key_sample, key_capture]:

            tracks_grouped = make_group_track(
                bigwigs,
                key,
                overlay=True,
                overlay_exclude=["subtraction", 
                                 "_vs_", 
                                 *re.split(r'[,;\s+]', P.PARAMS.get('compare_summary_methods', ['mean',]))],
            )

            trackdb.add_tracks(tracks_grouped.values())

        trackdb.validate()

        if P.PARAMS.get("hub_upload"):  # If the hub need to be uploaded to a server
            trackhub.upload.upload_hub(
                hub=hub, host=P.PARAMS["hub_url"], remote_dir=P.PARAMS["hub_dir"]
            )
        else:
            trackhub.upload.stage_hub(hub=hub, staging=P.PARAMS["hub_dir"])

    else:
        raise NotImplementedError("Custom genome not yet supported")



######################################
# Identify differential interactions #
######################################


@active_if(False)
@active_if(N_SAMPLES >= 4)
@follows(mkdir("ccanalyser_compare/differential"))
@transform(
    reporters_make_union_bedgraph,
    regex(r".*/(.*)\.raw\.tsv"),
    r"ccanalyser_compare/differential/\1.completed",
    extras=[r"\1"],
)
def identify_differential_interactions(infile, outfile, capture_name):

    if len(pd.read_csv(infile, sep="\t", nrows=5).columns) >= 4:

        output_prefix = outfile.replace(".log", "")

        statement = """ccanalyser
                       reporters 
                       differential
                    %(infile)s
                    -n %(capture_name)s
                    -c %(analysis_viewpoints)s
                    -o %(output_prefix)s
                    """

        P.run(
            statement,
            job_queue=P.PARAMS["pipeline_cluster_queue"],
            job_condaenv=P.PARAMS["conda_env"],
        )

    else:
        print("Not enough replicates for differential testing")
    
    touch_file(outfile)


##################
# Plot reporters #
##################


@active_if(VALID_PLOT_COORDINATES)
@follows(reporters_store_merged, mkdir("ccanalyser_analysis/heatmaps/"))
@transform(
    "ccanalyser_analysis/reporters/*.hdf5",
    regex(r"ccanalyser_analysis/reporters/(.*).hdf5"),
    r"ccanalyser_analysis/heatmaps/\1.completed",
)
def reporters_plot_heatmap(infile, outfile):
    """Plots a heatmap over a specified region"""

    if P.PARAMS.get("plot_normalisation"):
        norm = P.PARAMS["plot_normalisation"]
    else:
        norm_default = {
            "capture": "n_interactions",
            "tri": "n_rf_n_interactions",
            "tiled": "ice",
        }
        norm = norm_default[P.PARAMS["analysis_method"]]

    output_prefix = outfile.replace(".completed", "")

    resolutions = " -r ".join(re.split(r"[,;]\s*|\s+", str(P.PARAMS["plot_bin_size"])))

    statement = f"""ccanalyser reporters  plot
                    %(infile)s
                    -r %(resolutions)s
                    -c %(plot_coordinates)s
                    --normalisation %(norm)s
                    --cmap {P.PARAMS.get('plot_cmap', 'jet')}
                    --vmin {P.PARAMS.get('plot_min', '0')}
                    --vmax {P.PARAMS.get('plot_vmax', '1')}
                    -o %(output_prefix)s
                    > %(outfile)s"""

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    touch_file(outfile)


@follows(
    hub_make,
    reporters_plot_heatmap,
    reporters_make_union_bedgraph,
    identify_differential_interactions,
    reporters_make_comparison_bedgraph,
)
@originate(
    "pipeline_complete.txt",
)
def full(outfile):

    if os.path.exists("chrom_sizes.txt.tmp"):
        os.unlink("chrom_sizes.txt.tmp")

    if os.path.exists("ccanalyser_analysis/reporters/binners.pkl"):
        zap_file("ccanalyser_analysis/reporters/binners.pkl")

    touch_file(outfile)


if __name__ == "__main__":

    if (
        "-h" in sys.argv or "--help" in sys.argv
    ):  # If --help then just run the pipeline without setup
        P.main(sys.argv)
    elif not "make" in sys.argv:
        P.main(sys.argv)
    else:
        set_up_chromsizes()
        modify_pipeline_params_dict()
        check_user_supplied_paths()
        P.main(sys.argv)


# @transform(
#     merge_interactions,
#     regex(r".*/(.*).hdf5"),
#     r"ccanalyser_analysis/bedgraphs/\1.windowed.log",
#     extras=[r"\1"],
# )
# def make_bedgraph_windowed(infile, outfiles, sample_name):
#     """Extract reporters in bedgraph format from stored interactions"""

#     output_prefix = f"ccanalyser_analysis/bedgraphs/{sample_name}.windowed"

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

# @active_if(P.PARAMS.get("hub_url"))
# @follows(hub_make)
# @originate("hub_url.txt")
# def hub_write_path(outfile):
#     """Convinence task to write hub url to use for adding custom hub to UCSC genome browser"""

#     with open(outfile, "w") as w:
#         url = P.PARAMS["hub_url"].rstrip("/")
#         name_dir = P.PARAMS["hub_dir"].strip("/")
#         name_hubtxt = P.PARAMS["hub_name"] + ".hub.txt"

#         path_hubtxt = f"{url}/{name_dir}/{name_hubtxt}"

#         w.write(path_hubtxt)