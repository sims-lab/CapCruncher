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
from math import inf
import os
from posixpath import dirname
import re
import sys
import pickle
from cgatcore import pipeline as P
from cgatcore.iotools import touch_file, zap_file
import itertools
import warnings
import glob
import shutil
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
from capcruncher.tools.statistics import (
    collate_slice_data,
    collate_read_data,
    collate_cis_trans_data,
    collate_histogram_data,
    extract_trimming_stats,
)

from capcruncher.utils import is_on, is_none, is_valid_bed


##############################
#   Set-up global parameters #
##############################

# Set up global parameters dict
P.get_parameters("config.yml")


USE_BLACKLIST = is_valid_bed(
    P.PARAMS["analysis_optional_blacklist"]
)  # Determines if blacklist is used
VALID_PLOT_COORDINATES = is_valid_bed(
    P.PARAMS.get("plot_coordinates"), verbose=False
)  # Has valid plotting coordinate bed file
FASTQ_DEDUPLICATE = P.PARAMS.get(
    "deduplication_pre-dedup", False
)  # Turns on FASTQ deduplication
MAKE_HUB = is_on(P.PARAMS.get("hub_create"))  # Create a UCSC hub or not
N_SAMPLES = len(
    {re.match(r"(.*)_R*[12].fastq.*", fn).group(1) for fn in glob.glob("*.fastq*")}
)  # Determine the number of unique sample names


##############################
#  Pipeline set-up functions #
##############################


def check_config():
    """
    Checks that all essential configuration has been provided.
    """

    if not os.path.exists("config.yml"):
        raise OSError(
            "Configuration file: config.yml. Not present in working directory"
        )

    essential_keys = [
        "analysis_method",
        "analysis_viewpoints",
        "analysis_restriction_enzyme",
        "genome_name",
        "genome_fasta",
        "genome_aligner_index",
    ]

    for key in essential_keys:
        if not key in P.PARAMS:
            raise ValueError(
                f"No value provided for {key} in config.yml. Please correct this and re-run."
            )


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

    if P.PARAMS["genome_chrom_sizes"] and os.path.exists(
        P.PARAMS["genome_chrom_sizes"]
    ):
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


@follows(mkdir("capcruncher_preprocessing/restriction_enzyme_map/"))
@transform(
    P.PARAMS.get("genome_fasta"),
    regex(r".*/(.*).fa.*"),
    r"capcruncher_preprocessing/restriction_enzyme_map/genome.digest.bed.gz",
)
def genome_digest(infile, outfile):
    """
    In silco digestion of the genome to identify restriction fragment coordinates.

    Runs :ref:`capcruncher genome digest <CLI Documentation>`.

    """
    tmp = outfile.replace(".gz", "")
    statement_digest = " ".join(
        [
            "capcruncher",
            "genome",
            "digest",
            infile,
            "-r",
            P.PARAMS["analysis_restriction_enzyme"],
            "-o",
            tmp,
            "-l",
            f"{tmp}.log",
            "--sort",
        ]
    )

    statement_compress = " ".join(["pigz", "-p", "4", tmp])

    P.run(
        " && ".join([statement_digest, statement_compress]),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#############################
# Fastq file pre-processing #
#############################


@follows(mkdir("capcruncher_preprocessing"), mkdir("capcruncher_preprocessing/fastqc"))
@transform(
    "*.fastq*", regex(r"(.*).fastq.*"), r"capcruncher_preprocessing/fastqc/\1_fastqc.zip"
)
def fastq_qc(infile, outfile):
    """Runs fastqc on the input files to generate fastq statistics."""

    outdir = os.path.dirname(outfile)
    statement = " ".join(
        [
            "fastqc",
            infile,
            "-q",
            "-t",
            str(P.PARAMS["pipeline_n_cores"]),
            "--nogroup",
            "--outdir",
            outdir,
        ]
    )

    P.run(
        statement,
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("statistics"))
@merge(fastq_qc, "capcruncher_statistics/fastqc_report.html")
def fastq_multiqc(infile, outfile):
    """Collate fastqc reports into single report using multiqc"""

    basename = os.path.basename(outfile)
    dirname = os.path.dirname(outfile)

    statement_cleanup = " ".join(["rm", "-f", outfile])
    statement_export_1 = " ".join(["export", "LC_ALL=en_US.UTF-8"])
    statement_export_2 = " ".join(["export", "LANG=en_US.UTF-8"])
    statement_multiqc = " ".join(
        ["multiqc", "capcruncher_preprocessing/fastqc/", "-o", dirname, "-n", basename]
    )

    P.run(
        " && ".join(
            [
                statement_cleanup,
                statement_export_1,
                statement_export_2,
                statement_multiqc,
            ]
        ),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="2G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("capcruncher_preprocessing/split"))
@collate(
    "*.fastq.gz",
    regex(r"(.*)_R*[12].fastq.*"),
    r"capcruncher_preprocessing/split/\1.completed",
)
def fastq_split(infiles, outfile):
    """
    Splits the input fastq files into chunks for parallel processing

    Runs :ref:`capcruncher fastq split <CLI Documentation>`.

    """

    statement = [
        "capcruncher",
        "fastq",
        "split",
        " ".join(infiles),
        "-m",
        "unix",
        "-o",
        outfile.replace(".completed", ""),
        "-n",
        str(P.PARAMS.get("split_n_reads", 1e6)),
        "--no-gzip",
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Create sentinel file
    touch_file(outfile)


@active_if(FASTQ_DEDUPLICATE)
@follows(
    mkdir("capcruncher_preprocessing/deduplicated"),
    mkdir("capcruncher_preprocessing/deduplicated/deduplicated_ids"),
    fastq_split,
)
@collate(
    "capcruncher_preprocessing/split/*.fastq*",
    regex(r"capcruncher_preprocessing/split/(.*)_part(\d+)_[12].fastq(?:.gz)?"),
    r"capcruncher_preprocessing/deduplicated/deduplicated_ids/\1_\2.json.gz",
    extras=[r"\1", r"\2"],
)
def fastq_duplicates_parse(infiles, outfile, sample_name, part_no):

    """
    Parses fastq files into json format for sequence deduplication.

    Runs :ref:`capcruncher fastq deduplicate parse <CLI Documentation>`

    """

    statement = [
        "capcruncher",
        "fastq",
        "deduplicate",
        "parse",
        *[os.path.abspath(fn) for fn in infiles],
        "-o",
        outfile,
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="6G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@collate(
    fastq_duplicates_parse,
    regex(r"capcruncher_preprocessing/deduplicated/deduplicated_ids/(.*)_\d*.json.gz"),
    r"capcruncher_preprocessing/deduplicated/deduplicated_ids/\1.json.gz",
)
def fastq_duplicates_identify(infiles, outfile):

    """
    Identifies duplicate sequences from parsed fastq files in json format.

    Runs :ref:`capcruncher fastq deduplicate identify <CLI Documentation>`

    """

    statement = [
        "capcruncher",
        "fastq",
        "deduplicate",
        "identify",
        " ".join(infiles),
        "-o",
        outfile,
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory="32G",
        job_condaenv=P.PARAMS["conda_env"],
    )

    for fn in infiles:
        zap_file(fn)


@follows(
    fastq_duplicates_parse,
    fastq_duplicates_identify,
    mkdir("capcruncher_statistics/deduplication/data/"),
)
@collate(
    "capcruncher_preprocessing/split/*.fastq*",
    regex(r".*/(.*_part\d+)_[12].fastq(?:.gz)?"),
    r"capcruncher_preprocessing/deduplicated/\1.completed",
)
def fastq_duplicates_remove(infiles, outfile):

    """
    Removes duplicate read fragments identified from parsed fastq files.
    """

    sample = re.match(r".*/(.*)(_part\d+)_[12].fastq(?:.gz)?", infiles[0])
    sample_name = sample.group(1)
    sample_part = sample.group(2)
    output_prefix = outfile.replace(".completed", "")
    stats_prefix = (
        f"capcruncher_statistics/deduplication/data/{sample_name}{sample_part}"
    )

    if FASTQ_DEDUPLICATE:
        statement = " ".join(
            [
                "capcruncher",
                "fastq",
                "deduplicate",
                "remove",
                *infiles,
                "-d",
                f"capcruncher_preprocessing/deduplicated/deduplicated_ids/{sample_name}.json.gz",
                "--sample_name",
                sample_name,
                "--stats_prefix",
                stats_prefix,
                "-o",
                output_prefix,
            ]
        )

    else:
        statement = f"""ln -s $(pwd)/{infiles[0]} {output_prefix}_1.fastq &&
                        ln -s $(pwd)/{infiles[1]} {output_prefix}s_2.fastq &&
                        lc=$(cat {infiles[0]} | wc -l);
                        statsfile={stats_prefix}.deduplication.csv;
                        echo "stat,stat_type,read_type,read_number,stage,sample" > $statsfile;
                        echo -e $(($lc / 4)),reads_total,pe,0,deduplication,{sample_name} >> $statsfile;
                        echo -e $(($lc / 4)),reads_unique,pe,0,deduplication,{sample_name} >> $statsfile;
                        echo -e 0,reads_removed,pe,0,deduplication,{sample_name} >> $statsfile
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
    if FASTQ_DEDUPLICATE:
        for fn in infiles:
            zap_file(fn)


@follows(fastq_duplicates_remove)
@merge(
    "capcruncher_statistics/deduplication/data/*.csv",
    "capcruncher_statistics/deduplication/deduplication.reads.csv",
)
def stats_deduplication_collate(infiles, outfile):

    """Combines deduplication statistics from fastq file partitions."""

    stats_prefix = outfile.replace(".reads.csv", "")

    df_stats = collate_read_data(infiles)

    df_stats_read = df_stats.query('stat_type != "reads_removed"')

    df_stats.to_csv(f"{stats_prefix}.summary.csv", index=False)

    # Modified to enable more streamlined summary at final stage
    df_stats_read.to_csv(outfile, index=False)


@follows(
    mkdir("capcruncher_preprocessing/trimmed"),
    fastq_duplicates_remove,
    mkdir("capcruncher_statistics/trimming/data/"),
)
@collate(
    "capcruncher_preprocessing/deduplicated/*.fastq*",
    regex(r"capcruncher_preprocessing/deduplicated/(.*)_[12].fastq(?:.gz)?"),
    r"capcruncher_preprocessing/trimmed/\1.completed",
)
def fastq_trim(infiles, outfile):

    """Trim adaptor sequences from fastq files using trim_galore"""

    statement_trim = " ".join(
        [
            "trim_galore",
            " ".join(infiles),
            "-o",
            os.path.dirname(outfile),
            "--paired",
            "--cores",
            str(P.PARAMS.get("pipeline_n_cores", 1)),
            "--gzip",
            P.PARAMS.get("trim_options") or " ",
        ]
    )
    statement_stats_1 = " ".join(
        [
            "mv",
            f"capcruncher_preprocessing/trimmed/{os.path.basename(infiles[0])}_trimming_report.txt",
            "capcruncher_statistics/trimming/data",
        ]
    )
    statement_stats_2 = " ".join(
        [
            "mv",
            f"capcruncher_preprocessing/trimmed/{os.path.basename(infiles[1])}_trimming_report.txt",
            "capcruncher_statistics/trimming/data",
        ]
    )

    P.run(
        " && ".join([statement_trim, statement_stats_1, statement_stats_2]),
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
@merge(
    "capcruncher_statistics/trimming/data/*.txt",
    r"capcruncher_statistics/trimming/trimming.summary.csv",
)
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

    df_trimming_stats.to_csv(
        "capcruncher_statistics/trimming/trimming.summary.csv", index=False
    )


@follows(fastq_trim, mkdir("capcruncher_preprocessing/flashed"))
@collate(
    "capcruncher_preprocessing/trimmed/*.fq*",
    regex(r"capcruncher_preprocessing/trimmed/(.*)_[12]_.*.fq(?:.gz)?"),
    r"capcruncher_preprocessing/flashed/\1.completed",
)
def fastq_flash(infiles, outfile):

    """Combine overlapping paired-end reads using FLASh"""

    statement = [
        "flash",
        " ".join(infiles),
        "-o",
        outfile.replace(".completed", ""),
        "-t",
        str(P.PARAMS.get("pipeline_n_cores", 1)),
        "-z",
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    touch_file(outfile)


@follows(
    mkdir("capcruncher_preprocessing/digested"),
    fastq_flash,
    mkdir("capcruncher_statistics/digestion/data"),
)
@transform(
    "capcruncher_preprocessing/flashed/*.fastq.gz",
    regex(r"capcruncher_preprocessing/flashed/(.*).extendedFrags.fastq.gz"),
    r"capcruncher_preprocessing/digested/\1.flashed.fastq.gz",
)
def fastq_digest_combined(infile, outfile):

    """In silico restriction enzyme digest of combined (flashed) read pairs"""

    statement = [
        "capcruncher",
        "fastq",
        "digest",
        infile,
        "-o",
        outfile,
        "-m",
        "flashed",
        "-r",
        P.PARAMS["analysis_restriction_enzyme"],
        "--minimum_slice_length",
        "18",
        "--stats_prefix",
        f"capcruncher_statistics/digestion/data/{os.path.basename(outfile)}",
        "--sample_name",
        re.match(r".*/(.*?)_part.*", infile).group(1),
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=3,
        job_condaenv=P.PARAMS["conda_env"],
    )

    zap_file(infile)


@follows(fastq_flash, mkdir("capcruncher_statistics/digestion"))
@collate(
    "capcruncher_preprocessing/flashed/*.fastq.gz",
    regex(r"capcruncher_preprocessing/flashed/(.*).notCombined_[12].fastq.gz"),
    r"capcruncher_preprocessing/digested/\1.pe.fastq.gz",
)
def fastq_digest_non_combined(infiles, outfile):

    """In silico restriction enzyme digest of non-combined (non-flashed) read pairs"""

    statement = [
        "capcruncher",
        "fastq",
        "digest",
        " ".join(infiles),
        "-o",
        outfile,
        "-m",
        "pe",
        "-r",
        P.PARAMS["analysis_restriction_enzyme"],
        "--minimum_slice_length",
        "18",
        "--stats_prefix",
        f"capcruncher_statistics/digestion/data/{os.path.basename(outfile)}",
        "--sample_name",
        re.match(r".*/(.*?)_part.*", infiles[0]).group(1),
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=3,
        job_condaenv=P.PARAMS["conda_env"],
    )

    for fn in infiles:
        zap_file(fn)


@follows(fastq_digest_combined, fastq_digest_non_combined)
@merge(
    "capcruncher_statistics/digestion/data/*",
    "capcruncher_statistics/digestion/digestion.reads.csv",
)
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


@follows(mkdir("capcruncher_preprocessing"), fastq_preprocessing)
@transform(
    [fastq_digest_combined, fastq_digest_non_combined],
    regex(r"capcruncher_preprocessing/digested/(.*).fastq.gz"),
    r"capcruncher_preprocessing/aligned/\1.bam",
)
def fastq_alignment(infile, outfile):

    """Aligns in silico digested fastq files to the genome."""

    statement_align = " ".join(
        [
            P.PARAMS.get("align_aligner", "bowtie2"),
            P.PARAMS.get("align_index_flag") or " ",
            P.PARAMS["genome_aligner_index"],
            P.PARAMS.get("align_options") or "",
            infile,
        ]
    )

    statement_samtools_view = " ".join(["samtools", "view", "-b", "-S", ">", outfile])
    statement_samtools_sort = " ".join(
        [
            "samtools",
            "sort",
            outfile,
            "-o",
            f"{outfile}.sorted.bam",
            "-m",
            "2G",
            "-@",
            str(P.PARAMS["pipeline_n_cores"]),
        ]
    )
    statement_rename_bam = " ".join(["mv", "-f", f"{outfile}.sorted.bam", outfile])

    P.run(
        " ".join(
            [
                statement_align,
                "|",
                statement_samtools_view,
                "&&",
                statement_samtools_sort,
                "&&",
                statement_rename_bam,
            ]
        ),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory="4G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@collate(
    fastq_alignment,
    regex(r"capcruncher_preprocessing/aligned/(.*)_part\d+.*.bam"),
    r"capcruncher_preprocessing/aligned/\1.bam",
)
def alignments_merge(infiles, outfile):
    """
    Combines bam files (by flashed/non-flashed status and sample).

    This task simply provides an input for picard CollectAlignmentSummaryMetrics
    and is only used to provide overall mapping statistics. Fastq partitions
    are *not* combined at this stage.

    """

    statement = ["samtools", "merge", outfile, " ".join(infiles)]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@transform(alignments_merge, regex("(.*).bam"), r"\1.bam.bai")
def alignments_index(infile, outfile):

    """Indexes all bam files (both partitioned and merged)"""

    statement = [
        "samtools",
        "index",
        infile,
    ]

    P.run(
        " ".join(statement),
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


@originate("capcruncher_analysis/annotations/exclude.bed")
def annotate_make_exclusion_bed(outfile):

    """Generates exclusion window around each capture site"""

    if not is_valid_bed(P.PARAMS["analysis_viewpoints"]):
        raise ValueError("Viewpoints bed file is not a valid bed file")

    statement_bedtools_slop = " ".join(
        [
            "bedtools",
            "slop",
            "-i",
            P.PARAMS["analysis_viewpoints"],
            "-g",
            P.PARAMS["genome_chrom_sizes"],
            "-b",
            str(P.PARAMS["analysis_reporter_exclusion_zone"]),
        ]
    )
    statement_bedtools_subtract = " ".join(
        ["bedtools", "subtract", "-a", "-", "-b", P.PARAMS["analysis_viewpoints"]]
    )
    statement_sort = " ".join(["sort", "-k1,1", "-k2,2n", ">", outfile])

    P.run(
        "|".join(
            [statement_bedtools_slop, statement_bedtools_subtract, statement_sort]
        ),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@originate("capcruncher_analysis/annotations/viewpoints.bed")
def annotate_sort_viewpoints(outfile):

    """Sorts the capture oligos for bedtools intersect with --sorted option"""

    if not is_valid_bed(P.PARAMS["analysis_viewpoints"]):
        raise ValueError("Viewpoints bed file is not a valid bed file")

    statement = [
        "sort",
        "-k1,1",
        "-k2,2n",
        P.PARAMS["analysis_viewpoints"],
        ">",
        outfile,
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@originate("capcruncher_analysis/annotations/blacklist.bed")
def annotate_sort_blacklist(outfile):

    """Sorts the capture oligos for bedtools intersect with --sorted option"""

    if USE_BLACKLIST:
        statement = [
            "sort",
            "-k1,1",
            "-k2,2n",
            P.PARAMS["analysis_optional_blacklist"],
            ">",
            outfile,
        ]
    else:

        statement = ["touch", outfile]  # Make a blank file if no blacklist

    P.run(
        " ".join(statement),
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
                "fn": "capcruncher_preprocessing/restriction_enzyme_map/genome.digest.bed.gz",
                "action": "get",
                "fraction": 0.2,
            },
            {
                "name": "capture",
                "fn": "capcruncher_analysis/annotations/viewpoints.bed",
                "action": "get",
                "fraction": 0.9,
            },
            {
                "name": "exclusion",
                "fn": "capcruncher_analysis/annotations/exclude.bed",
                "action": "get",
                "fraction": 1e-9,
            },
            {
                "name": "exclusion_count",
                "fn": "capcruncher_analysis/annotations/exclude.bed",
                "action": "count",
                "fraction": 1e-9,
            },
            {
                "name": "capture_count",
                "fn": "capcruncher_analysis/annotations/viewpoints.bed",
                "action": "count",
                "fraction": 0.9,
            },
            {
                "name": "blacklist",
                "fn": "capcruncher_analysis/annotations/blacklist.bed",
                "action": "count",
                "fraction": 1e-9,
            },
        ]
    ),
    r"capcruncher_analysis/annotations/\1.annotations.tsv",
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

    flags = {"name": "-n", "fn": "-b", "action": "-a", "fraction": "-f"}
    statement_bamtobed = " ".join(["bedtools", "bamtobed", "-i", infile[0]])
    statement_sort = " ".join(["sort", "-k1,1", "-k2,2n"])
    statement_annotate = " ".join(
        [
            "capcruncher",
            "alignments",
            "annotate",
            "-",
            *[
                f"{flags[k]} {v}"
                for annotation in infile[1]
                for k, v in annotation.items()
            ],
            "-o",
            outfile,
            "--invalid_bed_action",
            "ignore",
            "-p",
            "1",
        ]
    )

    P.run(
        " ".join([statement_bamtobed, "|", statement_sort, "|", statement_annotate]),
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
    mkdir("capcruncher_analysis/reporters/unfiltered/"),
    mkdir("capcruncher_statistics/reporters/data"),
)
@transform(
    fastq_alignment,
    regex(r"capcruncher_preprocessing/aligned/(.*).bam"),
    add_inputs(r"capcruncher_analysis/annotations/\1.annotations.tsv"),
    r"capcruncher_analysis/reporters/unfiltered/\1.completed",
)
def alignments_filter(infiles, outfile):
    """Filteres slices and outputs reporter slices for each capture site"""

    bam, annotations = infiles
    sample = re.match(r".*/(.*)_(part\d+).(flashed|pe).bam", bam)
    sample_name = sample.group(1)
    sample_part = sample.group(2)
    sample_read_type = sample.group(3)

    output_prefix = outfile.replace(".completed", "")
    output_log_file = f"{output_prefix}.log"
    stats_prefix = f"capcruncher_statistics/reporters/data/{sample_name}_{sample_part}_{sample_read_type}"

    statement = [
        "capcruncher",
        "alignments",
        "filter",
        P.PARAMS["analysis_method"],
        "-b",
        bam,
        "-a",
        annotations,
        "-o",
        output_prefix,
        "--stats_prefix",
        stats_prefix,
        "--sample_name",
        sample_name,
        "--read_type",
        sample_read_type,
        ">",
        output_log_file,
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_memory=P.PARAMS["pipeline_memory"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Make sentinel file
    touch_file(outfile)

    # Zero annotations
    if not P.PARAMS.get("analysis_optional_keep_annotations", False):
        zap_file(annotations)


@follows(mkdir("capcruncher_analysis/reporters/collated"), alignments_filter)
@collate(
    "capcruncher_analysis/reporters/unfiltered/*.tsv",
    regex(
        r".*/(?P<sample>.*)_part\d+.(flashed|pe).(?P<capture>.*).(slices|fragments).tsv"
    ),
    r"capcruncher_analysis/reporters/collated/\1.\2.\3.\4.tsv",
    extras=[r"\1", r"\2", r"\3", r"\4"],
)
def reporters_collate(infiles, outfile, *grouping_args):

    """Concatenates identified reporters"""

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


@follows(alignments_filter, mkdir("capcruncher_analysis/reporters/deduplicated"))
@transform(
    reporters_collate,
    regex(r".*/(?P<sample>.*).(flashed|pe).(?P<capture>.*).fragments.tsv"),
    r"capcruncher_analysis/reporters/deduplicated/\1.\2.\3.json.gz",
    extras=[r"\2"],
)
def alignments_deduplicate_fragments(infile, outfile, read_type):

    """
    Identifies duplicate fragments with the same coordinates and order.
    """

    statement = [
        "capcruncher",
        "alignments",
        "deduplicate",
        "identify",
        infile,
        "--read_type",
        read_type,
        "-o",
        outfile,
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_memory="32G",
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(alignments_deduplicate_fragments)
@transform(
    reporters_collate,
    regex(
        r"capcruncher_analysis/reporters/collated/(.*)\.(flashed|pe)\.(.*)\.slices.tsv"
    ),
    add_inputs(r"capcruncher_analysis/reporters/deduplicated/\1.\2.\3.json.gz"),
    r"capcruncher_analysis/reporters/deduplicated/\1.\2.\3.slices.tsv",
    extras=[r"\1", r"\2", r"\3"],
)
def alignments_deduplicate_slices(
    infile, outfile, sample_name, read_type, capture_oligo
):

    """Removes reporters with duplicate coordinates"""

    slices, duplicated_ids = infile
    stats_prefix = f"capcruncher_statistics/reporters/data/{sample_name}_{read_type}_{capture_oligo}"

    statement = [
        "capcruncher",
        "alignments",
        "deduplicate",
        "remove",
        slices,
        "-d",
        duplicated_ids,
        "-o",
        outfile,
        "--stats_prefix",
        stats_prefix,
        "--sample_name",
        sample_name,
        "--read_type",
        read_type,
    ]

    P.run(
        " ".join(statement),
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
    r"capcruncher_analysis/reporters/\1.\2.tsv.gz",
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
@merge(
    "capcruncher_statistics/reporters/data/*",
    "capcruncher_statistics/reporters/reporters.reads.csv",
)
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
def post_capcruncher_analysis():
    """Reporters have been identified, deduplicated and collated by sample/capture probe"""


####################
# Reporter storage #
####################


@follows(mkdir("capcruncher_analysis/reporters/counts"))
@transform(
    alignments_deduplicate_collate,
    regex(r"capcruncher_analysis/reporters/(.*)\.(.*).tsv.gz"),
    r"capcruncher_analysis/reporters/counts/\1.\2.tsv.gz",
)
def reporters_count(infile, outfile):

    """Counts the number of interactions identified between reporter restriction fragments"""

    statement = [
        "capcruncher",
        "reporters",
        "count",
        infile,
        "-o",
        outfile,
        "--remove_exclusions",
        ">",
        f"{outfile}.log",
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=2,
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(mkdir("capcruncher_analysis/reporters/fragments"))
@transform(
    reporters_count,
    regex(r"capcruncher_analysis/reporters/counts/(.*)\.(.*)\.tsv.gz"),
    add_inputs(genome_digest),
    r"capcruncher_analysis/reporters/fragments/\1.\2.fragments.hdf5",
    extras=[r"\1", r"\2"],
)
def reporters_store_restriction_fragment(infile, outfile, sample_name, capture_name):

    """Stores restriction fragment interaction counts in cooler format"""

    counts, rf_map = infile
    output_prefix = outfile.replace(f".{capture_name}.fragments", "")

    statement = [
        "capcruncher",
        "reporters",
        "store",
        "fragments",
        counts,
        "-f",
        rf_map,
        "-g",
        P.PARAMS["genome_name"],
        "-n",
        capture_name,
        "-c",
        P.PARAMS["analysis_viewpoints"],
        "-o",
        output_prefix,
        "--suffix",
        "fragments",
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@follows(genome_digest, reporters_count)
@originate(r"capcruncher_analysis/reporters/binners.pkl")
def generate_bin_conversion_tables(outfile):
    """
    Converts restriction fragments to genomic bins.

    Binning restriction fragments into genomic bins takes a substantial
    amount of time and memory. To avoid repeatedly performing the same action,
    bin conversion tables are calculated once for each required resolution and
    then stored as a pickle file.

    """

    from capcruncher.tools.storage import GenomicBinner

    frags = pd.read_csv(
        "capcruncher_preprocessing/restriction_enzyme_map/genome.digest.bed.gz",
        sep="\t",
        names=["chrom", "start", "end", "name"],
    )

    binner_dict = dict()
    for bs in re.split(r"[,;]\s*|\s+", str(P.PARAMS["analysis_bin_size"])):
        gb = GenomicBinner(
            chromsizes=P.PARAMS["genome_chrom_sizes"], fragments=frags, binsize=int(bs)
        )
        bct = (
            gb.bin_conversion_table
        )  # Property is cached so need to call it to make sure it is present.
        binner_dict[int(bs)] = gb

    with open("capcruncher_analysis/reporters/binners.pkl", "wb") as w:
        pickle.dump(binner_dict, w)


@active_if(P.PARAMS.get("analysis_bin_size"))
@follows(generate_bin_conversion_tables, mkdir("capcruncher_analysis/reporters/binned/"))
@transform(
    reporters_store_restriction_fragment,
    regex(r"capcruncher_analysis/reporters/fragments/(.*)\.(.*)\.fragments\.hdf5"),
    add_inputs(generate_bin_conversion_tables),
    r"capcruncher_analysis/reporters/binned/\1.\2.completed",
    extras=[r"\2"],
)
def reporters_store_binned(infile, outfile, capture_name):

    """
    Converts a cooler file of restriction fragments to even genomic bins.
    """

    clr, conversion_tables = infile
    statement = [
        "capcruncher",
        "reporters",
        "store",
        "bins",
        clr,
        *[
            f"-b {bin_size}"
            for bin_size in re.split(r"[,;]\s*|\s+", str(P.PARAMS["analysis_bin_size"]))
        ],
        "--conversion_tables",
        conversion_tables,
        "--normalise",
        "-p",
        str(P.PARAMS["pipeline_n_cores"]),
        "-o",
        outfile.replace(f".{capture_name}.completed", ""),
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_threads=P.PARAMS["pipeline_n_cores"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    # Make sentinel file
    touch_file(outfile)


@follows(reporters_store_restriction_fragment, reporters_store_binned)
@collate(
    [
        "capcruncher_analysis/reporters/fragments/*.hdf5",
        "capcruncher_analysis/reporters/binned/*.hdf5",
    ],
    regex(r".*/(.*)\.(.*)\.(?:fragments|\d+)\.hdf5"),
    r"capcruncher_analysis/reporters/\1.hdf5",
    extras=[r"\1"],
)
def reporters_store_merged(infiles, outfile, sample_name):

    """Combines cooler files together"""

    statement = [
        "capcruncher",
        "reporters",
        "store",
        "merge",
        " ".join(infiles),
        "-o",
        outfile,
    ]

    P.run(
        " ".join(statement),
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
        stats_trim_collate,
        stats_digestion_collate,
        stats_alignment_filtering_collate,
    ],
    "capcruncher_statistics/run_statistics.csv",
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
    "capcruncher_statistics/visualise_statistics.html",
)
def pipeline_make_report(infile, outfile):
    """Run jupyter notebook for reporting and plotting pipeline statistics"""

    path_script = __file__
    path_script_dir = os.path.dirname(path_script)
    path_nb_dir = os.path.dirname(path_script_dir)

    statement_clean = " ".join(["rm", outfile.replace(".html", "*"), "-f"])

    statement_papermill = " ".join(
        [
            "papermill",
            "-k",
            "python3",
            "-p",
            "directory",
            "$(pwd)/capcruncher_statistics/",
            f"{path_nb_dir}/visualise_statistics.ipynb",
            outfile.replace(".html", ".ipynb"),
        ]
    )

    statement_nbconvert = " ".join(
        [
            "jupyter",
            "nbconvert",
            "--no-input",
            "--to html",
            outfile.replace(".html", ".ipynb"),
            outfile,
        ]
    )

    P.run(
        " && ".join([statement_clean, statement_papermill, statement_nbconvert]),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#####################
# Reporter pileups  #
#####################


@follows(mkdir("capcruncher_analysis/bedgraphs"))
@transform(
    reporters_store_merged,
    regex(r".*/(.*).hdf5"),
    r"capcruncher_analysis/bedgraphs/\1.raw.completed",
    extras=[r"\1"],
)
def reporters_make_bedgraph(infile, outfile, sample_name):
    """Extract reporters in bedgraph format from stored interactions"""

    output_prefix = f"capcruncher_analysis/bedgraphs/{sample_name}.raw"

    statement = ["capcruncher", "reporters", "pileup", infile, "-o", output_prefix]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    touch_file(outfile)


@transform(
    reporters_store_merged,
    regex(r".*/(.*).hdf5"),
    r"capcruncher_analysis/bedgraphs/\1.normalised.completed",
    extras=[r"\1"],
)
def reporters_make_bedgraph_normalised(infile, outfile, sample_name):
    """
    Extract reporters in bedgraph format from stored interactions.

    In addition to generating a bedgraph this task also normalises the counts
    by the number of cis interactions identified to enable cross sample comparisons.

    """

    output_prefix = f"capcruncher_analysis/bedgraphs/{sample_name}.normalised"

    statement = [
        "capcruncher",
        "reporters",
        "pileup",
        infile,
        "-o",
        output_prefix,
        "--normalise",
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    touch_file(outfile)


@active_if(N_SAMPLES >= 2)
@follows(
    mkdir("capcruncher_compare/bedgraphs_union"),
    reporters_make_bedgraph,
    reporters_make_bedgraph_normalised,
)
@collate(
    "capcruncher_analysis/bedgraphs/*.bedgraph",
    regex(r".*/(?:.*)\.(raw|normalised|windowed)\.(.*).bedgraph"),
    r"capcruncher_compare/bedgraphs_union/\2.\1.tsv",
    extras=[r"\1", r"\2"],
)
def reporters_make_union_bedgraph(infiles, outfile, normalisation_type, capture_name):

    """
    Collates bedgraphs by capture probe into a single file for comparison.

    See `bedtools unionbedg <https://bedtools.readthedocs.io/en/latest/content/tools/unionbedg.html>`_
    for more details.

    """
    sample_names = [
        re.match(r".*/(.*)\.(.*)\.(?:.*).bedgraph(?:.gz)?", fn).group(1)
        for fn in infiles
    ]

    statement = [
        "bedtools",
        "unionbedg",
        "-i",
        " ".join(infiles),
        "-header",
        "-names",
        " ".join(sample_names),
        ">",
        outfile,
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(N_SAMPLES >= 2)
@follows(
    mkdir("capcruncher_compare/bedgraphs_comparison/"), reporters_make_union_bedgraph
)
@transform(
    reporters_make_union_bedgraph,
    regex(r"capcruncher_compare/bedgraphs_union/(.*)\.normalised\.tsv"),
    r"capcruncher_compare/bedgraphs_comparison/\1.completed",
    extras=[r"\1"],
)
def reporters_make_comparison_bedgraph(infile, outfile, viewpoint):

    import numpy as np

    df_bdg = pd.read_csv(infile, sep="\t")
    dir_output = os.path.dirname(outfile)

    summary_methods = [
        m
        for m in re.split(r"[,;\s+]", P.PARAMS.get("compare_summary_methods", "mean,"))
        if m
    ]
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
            a_summary = pd.Series(
                df_a.pipe(summary_functions[summary_method], axis=1),
                name=summary_method,
            )
            b_summary = pd.Series(
                df_b.pipe(summary_functions[summary_method], axis=1),
                name=summary_method,
            )

            df_a_bdg = pd.concat([df_bdg.iloc[:, :3], a_summary], axis=1)
            df_b_bdg = pd.concat([df_bdg.iloc[:, :3], b_summary], axis=1)
            df_subtraction_bdg = pd.concat(
                [df_bdg.iloc[:, :3], a_summary - b_summary], axis=1
            )

            df_a_bdg.to_csv(
                f"{dir_output}/{a}.{summary_method}-summary.{viewpoint}.bedgraph",
                sep="\t",
                header=False,
                index=None,
            )

            df_b_bdg.to_csv(
                f"{dir_output}/{b}.{summary_method}-summary.{viewpoint}.bedgraph",
                sep="\t",
                header=False,
                index=None,
            )

            df_subtraction_bdg.to_csv(
                f"{dir_output}/{a}_vs_{b}.{summary_method}-subtraction.{viewpoint}.bedgraph",
                sep="\t",
                index=None,
                header=False,
            )

    touch_file(outfile)


@follows(
    mkdir("capcruncher_analysis/bigwigs"),
    reporters_make_bedgraph,
    reporters_make_bedgraph_normalised,
    reporters_make_comparison_bedgraph,
)
@transform(
    [
        "capcruncher_analysis/bedgraphs/*",
        "capcruncher_compare/bedgraphs_comparison/*.bedgraph",
    ],
    regex(r".*/(.*).bedgraph"),
    r"capcruncher_analysis/bigwigs/\1.bigWig",
)
def reporters_make_bigwig(infile, outfile):
    """Uses UCSC tools bedGraphToBigWig to generate bigWigs for each bedgraph"""

    tmp = f"{outfile}.tmp"
    statement_sort = " ".join(["sort", "-k1,1", "-k2,2n", infile, ">", tmp])
    statement_bdgtobw = " ".join(
        ["bedGraphToBigWig", tmp, P.PARAMS["genome_chrom_sizes"], outfile]
    )
    statement_cleanup = " ".join(["rm", tmp])

    P.run(
        " && ".join([statement_sort, statement_bdgtobw, statement_cleanup]),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


#######################
# UCSC hub generation #
#######################


@mkdir("capcruncher_analysis/viewpoints/")
@transform(
    annotate_sort_viewpoints,
    regex(r".*/(.*).bed"),
    r"capcruncher_analysis/viewpoints/\1.bigBed",
)
def viewpoints_to_bigbed(infile, outfile):

    statement = ["bedToBigBed", infile, P.PARAMS["genome_chrom_sizes"], outfile]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )


@active_if(MAKE_HUB)
@merge(
    [reporters_make_bigwig, viewpoints_to_bigbed, pipeline_make_report],
    os.path.join(
        P.PARAMS.get("hub_dir", ""), P.PARAMS.get("hub_name", "") + ".hub.txt"
    ),
)
def hub_make(infiles, outfile):
    """Creates a ucsc hub from the pipeline output"""

    import trackhub
    import seaborn as sns
    from capcruncher.utils import categorise_tracks

    # Extract statistics
    stats_report = [fn for fn in infiles if fn.endswith(".html")][0]

    # Create a dataframe with bigwig attributes and paths
    bigwigs = [fn for fn in infiles if fn.endswith(".bigWig")]
    df_bigwigs = (
        pd.Series(bigwigs)
        .to_frame("fn")
        .assign(basename=lambda df: df["fn"].apply(os.path.basename))
    )
    attributes = df_bigwigs["basename"].str.extract(
        r"(?P<samplename>.*?)\.(?P<method>.*?)\.(?P<viewpoint>.*?)\.(?P<filetype>.*)"
    )
    df_bigwigs = (
        df_bigwigs.join(attributes)
        .assign(track_categories=lambda df: categorise_tracks(df["method"]))
        .sort_values(["samplename", "method", "viewpoint"])
    )

    # Create a hub
    hub = trackhub.Hub(
        hub=P.PARAMS["hub_name"],
        short_label=P.PARAMS.get("hub_short", P.PARAMS["hub_name"]),
        long_label=P.PARAMS.get("hub_long", P.PARAMS["hub_name"]),
        email=P.PARAMS["hub_email"],
    )

    ## Need to make an assembly hub if this is a custom genome
    if P.PARAMS.get("genome_custom"):
        genome = trackhub.Assembly(
            genome=P.PARAMS["genome_name"],
            twobit_file=P.PARAMS["genome_twobit"],
            organism=P.PARAMS["genome_organism"],
            defaultPos=P.PARAMS.get("hub_default_position", "chr1:1000-2000"),
        )
        groups_file = trackhub.GroupsFile(
            [
                trackhub.GroupDefinition(
                    name=P.PARAMS["hub_name"], priority=1, default_is_closed=False
                ),
            ]
        )
        genome.add_groups(groups_file)

    else:
        genome = trackhub.Genome(P.PARAMS["genome_name"])
        groups_file = None

    # Create genomes file
    genomes_file = trackhub.GenomesFile()

    # Create trackdb
    trackdb = trackhub.TrackDb()

    # Add these to the hub
    hub.add_genomes_file(genomes_file)
    genome.add_trackdb(trackdb)
    genomes_file.add_genome(genome)

    # Extract groups for generating composite tracks
    unique_samples = df_bigwigs["samplename"].unique()
    unique_viewpoints = df_bigwigs["viewpoint"].unique()
    unique_comparison_methods = df_bigwigs["method"].unique()

    subgroup_vp = trackhub.SubGroupDefinition(
        name="viewpoint",
        label="Viewpoint",
        mapping={n.lower(): n for n in unique_viewpoints},
    )
    subgroup_sample = trackhub.SubGroupDefinition(
        name="samplename",
        label="Sample_Name",
        mapping={n.lower(): n for n in unique_samples},
    )
    subgroup_method = trackhub.SubGroupDefinition(
        name="summary_method",
        label="Summary_Method",
        mapping={
            n.split("-")[0]: n.split("-")[0]
            for n in unique_comparison_methods
        },
    )

    # Generate a color mapping based on sample names
    colors = sns.color_palette("hls", len(unique_samples))
    color_mapping = dict(zip(unique_samples, colors))

    #####################
    # Add tracks to hub #
    #####################

    for category_name, df in df_bigwigs.groupby("track_categories"):

        composite = trackhub.CompositeTrack(
            name=category_name,
            short_label=category_name,
            dimensions="dimX=samplename dimY=viewpoint dimA=summary_method",
            sortOrder="samplename=+ viewpoint=+ summary_method=+",
            tracktype="bigWig",
            visibility="hide",
            dragAndDrop="subTracks",
            allButtonPair="off",
        )

        # Only add a group if this is an assembly hub
        if groups_file:
            composite.add_params(group=P.PARAMS["hub_name"])

        composite.add_subgroups([subgroup_vp, subgroup_sample, subgroup_method])
        # composite.add_params(html=os.path.basename(stats_report))

        for bw in df.itertuples():
            t = trackhub.Track(
                name=f'{bw.samplename}_{bw.viewpoint}_{bw.method.replace("-summary", "")}',
                source=bw.fn,
                autoScale="off",
                tracktype="bigWig",
                windowingFunction="maximum",
                subgroups={
                    "viewpoint": bw.viewpoint.lower(),
                    "samplename": bw.samplename.lower(),
                    "summary_method": bw.method.split("-")[0],
                },
                color=",".join(
                    [str(int(x * 255)) for x in color_mapping[bw.samplename]]
                ),
            )

            # Only add a group if this is an assembly hub
            if groups_file:
                t.add_params(group=P.PARAMS["hub_name"])

            composite.add_subtrack(t)

        trackdb.add_tracks(composite)

    # Add viewpoints to hub
    for bb in [fn for fn in infiles if fn.endswith(".bigBed")]:

        t = trackhub.Track(
            name=os.path.basename(bb).replace(".bigBed", "").capitalize(),
            source=bb,
            tracktype="bigBed",
        )

        if genomes_file:
            t.add_params(group=P.PARAMS["hub_name"])

        trackdb.add_tracks(t)

    #############
    # Stage hub #
    #############

    staging_tmp_dir = "hub_tmp_dir"
    
    # Stage the hub
    trackhub.upload.stage_hub(hub=hub, staging=staging_tmp_dir)

    # Edit the hub.txt file to include the stats report as descriptionUrl
    with open(os.path.join(staging_tmp_dir, f'{P.PARAMS["hub_name"]}.hub.txt'), 'a') as hubtxt:
        hubtxt.write('\n')
        hubtxt.write(f'descriptionUrl {P.PARAMS["genome_name"]}/{os.path.basename(stats_report)}\n')

    # Copy to the new location
    shutil.copytree(
        "hub_tmp_dir",
        P.PARAMS["hub_dir"],
        dirs_exist_ok=True,
        symlinks=P.PARAMS.get("hub_symlink", False),
    )

    # Delete the staged hub
    shutil.rmtree("hub_tmp_dir")

    # Copy the stats report to the correct location
    shutil.copy(
        stats_report,
        os.path.join(
            P.PARAMS["hub_dir"], P.PARAMS["genome_name"], os.path.basename(stats_report)
        ),
    )



######################################
# Identify differential interactions #
######################################


@active_if(False)
@active_if(N_SAMPLES >= 4)
@follows(mkdir("capcruncher_compare/differential"))
@transform(
    reporters_make_union_bedgraph,
    regex(r".*/(.*)\.raw\.tsv"),
    r"capcruncher_compare/differential/\1.completed",
    extras=[r"\1"],
)
def identify_differential_interactions(infile, outfile, capture_name):

    if len(pd.read_csv(infile, sep="\t", nrows=5).columns) >= 4:

        output_prefix = outfile.replace(".log", "")

        statement = [
            "capcruncher",
            "reporters",
            "differential",
            infile,
            "-n",
            capture_name,
            "-c",
            P.PARAMS["analysis_viewpoints"],
            "-o",
            outfile.replace(".log", ""),
        ]

        P.run(
            " ".join(statement),
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
@follows(reporters_store_merged, mkdir("capcruncher_analysis/heatmaps/"))
@transform(
    "capcruncher_analysis/reporters/*.hdf5",
    regex(r"capcruncher_analysis/reporters/(.*).hdf5"),
    r"capcruncher_analysis/heatmaps/\1.completed",
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

    statement = [
        "capcruncher",
        "reporters",
        "plot",
        infile,
        "-r",
        resolutions,
        "-c",
        P.PARAMS["plot_coordinates"],
        "--normalisation",
        norm,
        "--cmap",
        P.PARAMS.get("plot_cmap", "jet"),
        "--vmin",
        P.PARAMS.get("plot_min", "0"),
        "--vmax",
        P.PARAMS.get("plot_vmax", "1"),
        "-o",
        output_prefix,
    ]

    P.run(
        " ".join(statement),
        job_queue=P.PARAMS["pipeline_cluster_queue"],
        job_condaenv=P.PARAMS["conda_env"],
    )

    touch_file(outfile)


@follows(
    pipeline_make_report,
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

    if os.path.exists("capcruncher_analysis/reporters/binners.pkl"):
        zap_file("capcruncher_analysis/reporters/binners.pkl")

    touch_file(outfile)


if __name__ == "__main__":

    if (
        "-h" in sys.argv or "--help" in sys.argv
    ):  # If --help then just run the pipeline without setup
        P.main(sys.argv)
    elif not "make" in sys.argv:
        P.main(sys.argv)
    else:
        check_config()
        set_up_chromsizes()
        modify_pipeline_params_dict()
        check_user_supplied_paths()
        P.main(sys.argv)
