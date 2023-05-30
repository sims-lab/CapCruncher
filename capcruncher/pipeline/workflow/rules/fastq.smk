import os
import pathlib
import json
import re


def get_parts(wc, sample_name=None):

    if not sample_name:
        sample_name = wc.sample

    checkpoint_output = checkpoints.split.get(sample=sample_name).output[0]

    parts = set()
    pattern = re.compile(r"(.*?)_part(\d+)_[12].fastq.gz")
    for fn in pathlib.Path(checkpoint_output).glob("*.fastq.gz"):
        parts.add(pattern.match(fn.name).group(2))

    return list(parts)


def get_fastq_parts(wc):
    return expand(
        "capcruncher_output/fastq/split/{{sample}}/{{sample}}_part{part}_{read}.fastq.gz",
        part=get_parts(wc),
        read=["1", "2"],
    )


def get_pickles(wc):
    return expand(
        "capcruncher_output/fastq/deduplicated/{{sample}}/{{sample}}_{part}.pkl",
        part=get_parts(wc),
    )


checkpoint split:
    input:
        fq1=lambda wc: FASTQ_SAMPLES.translation[f"{wc.sample}_1.fastq.gz"],
        fq2=lambda wc: FASTQ_SAMPLES.translation[f"{wc.sample}_2.fastq.gz"],
    output:
        directory("capcruncher_output/fastq/split/{sample}"),
    threads: config["pipeline"].get("n_cores", 1)
    params:
        prefix="capcruncher_output/fastq/split/{sample}/{sample}",
        n_reads=str(config["split"].get("n_reads", 1e6)),
    log:
        "capcruncher_output/logs/split/{sample}.log",
    shell:
        """
        mkdir {output} && \
        capcruncher \
        fastq \
        split \
        {input.fq1} \
        {input.fq2} \
        -m \
        unix \
        -o \
        {params.prefix} \
        -n \
        {params.n_reads} \
        --gzip \
        -p \
        {threads} \
        > {log} 2>&1
        """


rule deduplication_parse:
    input:
        fq1="capcruncher_output/fastq/split/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_output/fastq/split/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        temp("capcruncher_output/fastq/deduplicated/{sample}/{sample}_{part}.pkl"),
    log:
        "capcruncher_output/logs/deduplication_fastq/parse/{sample}_part{part}.log",
    shell:
        """
        capcruncher \
        fastq \
        deduplicate \
        parse \
        {input.fq1} \
        {input.fq2} \
        -o \
        {output} \
        > {log} 2>&1
        """


rule deduplication_identify:
    input:
        hashed_reads=get_pickles,
    output:
        temp("capcruncher_output/fastq/deduplicated/{sample}/{sample}.pkl"),
    log:
        "capcruncher_output/logs/deduplication_fastq/identify/{sample}.log",
    threads: 3
    shell:
        """
        capcruncher \
        fastq \
        deduplicate \
        identify \
        {input.hashed_reads} \
        -o \
        {output}
        > {log} 2>&1
        """


rule deduplication_remove:
    input:
        fq1="capcruncher_output/fastq/split/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_output/fastq/split/{sample}/{sample}_part{part}_2.fastq.gz",
        ids_duplicated="capcruncher_output/fastq/deduplicated/{sample}/{sample}.pkl",
    output:
        fq1=temp(
            "capcruncher_output/fastq/deduplicated/{sample}/{sample}_part{part}_1.fastq.gz"
        ),
        fq2=temp(
            "capcruncher_output/fastq/deduplicated/{sample}/{sample}_part{part}_2.fastq.gz"
        ),
        stats="capcruncher_output/statistics/deduplication/data/{sample}_part{part}.deduplication.csv",
    params:
        prefix_fastq="capcruncher_output/fastq/deduplicated/{sample}/{sample}_part{part}",
        prefix_stats="capcruncher_output/statistics/deduplication/data/{sample}_part{part}",
    log:
        "capcruncher_output/logs/deduplication_fastq/remove/{sample}_part{part}.log",
    threads: 4
    shell:
        """
        capcruncher \
        fastq \
        deduplicate \
        remove \
        {input.fq1} \
        {input.fq2} \
        -d \
        {input.ids_duplicated} \
        --sample-name \
        {wildcards.sample} \
        --stats-prefix \
        {params.prefix_stats} \
        -o \
        {params.prefix_fastq} \
        -p \
        {threads} \
        --hash-read-name \
        --gzip \
        > {log} 2>&1
        """


rule trim:
    input:
        fq1=rules.deduplication_remove.output.fq1,
        fq2=rules.deduplication_remove.output.fq2,
    output:
        trimmed1=temp(
            "capcruncher_output/fastq/trimmed/{sample}/{sample}_part{part}_1.fastq.gz"
        ),
        trimmed2=temp(
            "capcruncher_output/fastq/trimmed/{sample}/{sample}_part{part}_2.fastq.gz"
        ),
    params:
        outdir="capcruncher_output/fastq/trimmed/{sample}/",
    threads: 4
    log:
        "capcruncher_output/logs/trimming/{sample}_{part}.log",
    shell:
        """
           trim_galore --cores {threads} --trim-n --paired --output_dir {params.outdir} {input.fq1} {input.fq2} >> {log} 2>&1 &&
           mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_1_val_1.fq.gz {output.trimmed1} &&
           mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_2_val_2.fq.gz {output.trimmed2}
        """


rule flash:
    input:
        fq1=rules.trim.output.trimmed1,
        fq2=rules.trim.output.trimmed2,
    output:
        flashed=temp(
            "capcruncher_output/fastq/flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz"
        ),
        pe1=temp(
            "capcruncher_output/fastq/flashed/{sample}/{sample}_part{part}.notCombined_1.fastq.gz"
        ),
        pe2=temp(
            "capcruncher_output/fastq/flashed/{sample}/{sample}_part{part}.notCombined_2.fastq.gz"
        ),
    params:
        outdir="capcruncher_output/fastq/flashed/{sample}/{sample}_part{part}",
    threads: 12
    log:
        "capcruncher_output/logs/flash/{sample}_{part}.log",
    shell:
        """
        flash {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z --compress-prog-args pigz > {log} 2>&1
        """


rule digest_flashed_combined:
    input:
        flashed=rules.flash.output.flashed,
    output:
        digested=temp(
            "capcruncher_output/fastq/digested/{sample}/{sample}_part{part}_flashed.fastq.gz"
        ),
        stats_read="capcruncher_output/statistics/digestion/data/{sample}_part{part}_flashed.digestion.read.summary.csv",
        stats_unfiltered="capcruncher_output/statistics/digestion/data/{sample}_part{part}_flashed.digestion.unfiltered.histogram.csv",
        stats_filtered="capcruncher_output/statistics/digestion/data/{sample}_part{part}_flashed.digestion.filtered.histogram.csv",
    params:
        prefix_stats="capcruncher_output/statistics/digestion/data/{sample}_part{part}_flashed",
        restriction_site=config["analysis"]["restriction_enzyme"],
    threads: 8
    log:
        "capcruncher_output/logs/digestion/{sample}_{part}.log",
    shell:
        """
        capcruncher \
        fastq \
        digest \
        {input.flashed} \
        -o \
        {output.digested} \
        -p \
        {threads} \
        -m \
        flashed \
        -r \
        {params.restriction_site} \
        --minimum_slice_length \
        18 \
        --stats-prefix \
        {params.prefix_stats} \
        --sample-name \
        {wildcards.sample} \
        > {log} 2>&1
        """


rule digest_flashed_pe:
    input:
        pe1=rules.flash.output.pe1,
        pe2=rules.flash.output.pe2,
    output:
        digested=temp(
            "capcruncher_output/fastq/digested/{sample}/{sample}_part{part}_pe.fastq.gz"
        ),
        stats_read="capcruncher_output/statistics/digestion/data/{sample}_part{part}_pe.digestion.read.summary.csv",
        stats_unfiltered="capcruncher_output/statistics/digestion/data/{sample}_part{part}_pe.digestion.unfiltered.histogram.csv",
        stats_filtered="capcruncher_output/statistics/digestion/data/{sample}_part{part}_pe.digestion.filtered.histogram.csv",
    params:
        prefix_stats="capcruncher_output/statistics/digestion/data/{sample}_part{part}_pe",
        restriction_site=config["analysis"]["restriction_enzyme"],
    threads: 8
    log:
        "capcruncher_output/logs/digestion/{sample}_{part}.log",
    shell:
        """
        capcruncher \
        fastq \
        digest \
        {input.pe1} \
        {input.pe2} \
        -o \
        {output.digested} \
        -p \
        {threads} \
        -m \
        pe \
        -r \
        {params.restriction_site} \
        --minimum_slice_length \
        18 \
        --stats-prefix \
        {params.prefix_stats} \
        --sample-name \
        {wildcards.sample} \
        > {log} 2>&1
        """
