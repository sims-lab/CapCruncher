import os
import pathlib
import json
import re
from typing import Literal
import capcruncher.pipeline.utils


def get_split_1_parts(wildcards):

    import pathlib
    import json
    import re

    outdir = checkpoints.split.get(**wildcards).output[0]
    fq_files = pathlib.Path(outdir).glob("*.fastq.gz")
    parts = sorted(
        set(
            [
                int(re.search(r"part(\d+)", f.name).group(1))
                for f in fq_files
                if re.search(r"part(\d+)", f.name)
            ]
        )
    )

    return parts


def get_pickles(wc):
    return expand(
        "capcruncher_output/interim/fastq/deduplicated/{{sample}}/{{sample}}_{part}.pkl",
        part=get_split_1_parts(wc),
    )


def get_fastq_split_1(wildcards):
    return {
        f"fq{read}": expand(
            "capcruncher_output/interim/fastq/split/{{sample}}/{{sample}}_part{part}_{read}.fastq.gz",
            part=get_split_1_parts(wildcards),
            read=[read],
        )
        for read in ["1", "2"]
    }


def get_deduplicated_fastq_pair(wildcards):
    import pathlib

    input_dir = checkpoints.deduplication.get(**wildcards).output[0]

    fq = {
        f"fq{read}": f"{input_dir.rstrip('/')}/{wildcards.sample}_part{wildcards.part}_{read}.fastq.gz"
        for read in ["1", "2"]
    }

    if pathlib.Path(fq["fq1"]).exists() and pathlib.Path(fq["fq2"]).exists():
        return fq
    else:
        return {"fq1": [], "fq2": []}


def get_flashed_fastq(wildcards):
    import pathlib

    fq = [
        f"capcruncher_output/interim/fastq/flashed/{wildcards.sample}/{wildcards.sample}_part{part}.extendedFrags.fastq.gz"
        for part in get_split_1_parts(wildcards)
    ]
    return fq


def get_pe_fastq(wildcards):
    fq = [
        f"capcruncher_output/interim/fastq/flashed/{wildcards.sample}/{wildcards.sample}_part{part}.notCombined_{read}.fastq.gz"
        for part in get_split_1_parts(wildcards)
        for read in ["1", "2"]
    ]
    return fq


def get_rebalanced_parts(wc, combined: Literal["flashed", "pe"], sample: str = None):
    if not sample:
        sample = wc.sample

    if combined == "flashed":
        checkpoint_output = checkpoints.rebalance_partitions_combined.get(
            sample=sample
        ).output[0]
        parts = glob_wildcards(
            "capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/{sample_name}_part{part}_flashed_1.fastq.gz"
        ).part
    elif combined == "pe":
        checkpoint_output = checkpoints.rebalance_partitions_pe.get(
            sample=sample
        ).output[0]
        parts = glob_wildcards(
            "capcruncher_output/interim/fastq/rebalanced/{sample}/pe/{sample_name}_part{part}_pe_1.fastq.gz"
        ).part

    else:
        raise ValueError(f"Unknown combined type {combined}")

    return set(parts)


def get_rebalanced_fastq_combined(wc):
    checkpoint_output = checkpoints.rebalance_partitions_combined.get(**wc).output[0]
    return f"capcruncher_output/interim/fastq/rebalanced/{wc.sample}/flashed/{wc.sample}_part{wc.part}_flashed_1.fastq.gz"


def get_rebalanced_fastq_pe(wc):
    checkpoint_output = checkpoints.rebalance_partitions_pe.get(
        **wc,
    ).output[0]
    return {
        "pe1": f"capcruncher_output/interim/fastq/rebalanced/{wc.sample}/pe/{wc.sample}_part{wc.part}_pe_1.fastq.gz",
        "pe2": f"capcruncher_output/interim/fastq/rebalanced/{wc.sample}/pe/{wc.sample}_part{wc.part}_pe_2.fastq.gz",
    }


def get_deduplicated_fastq(wc):
    checkpoint_output = checkpoints.deduplication.get(sample=wc.sample).output[0]
    return {
        "fq1": f"capcruncher_output/interim/fastq/deduplicated/{wc.sample}/{wc.sample}_part{wc.part}_1.fastq.gz",
        "fq2": f"capcruncher_output/interim/fastq/deduplicated/{wc.sample}/{wc.sample}_part{wc.part}_2.fastq.gz",
    }


def separate_pe_fastq(wc):
    return {
        1: expand(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_{read}.fastq.gz",
            sample=wc.sample,
            part=get_split_1_parts(wc),
            read=["1"],
        ),
        2: expand(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_{read}.fastq.gz",
            sample=wc.sample,
            part=get_split_1_parts(wc),
            read=["2"],
        ),
    }


rule fastq_rename:
    input:
        fq1=lambda wc: FASTQ_SAMPLES.translation[f"{wc.sample}_1.fastq.gz"],
        fq2=lambda wc: FASTQ_SAMPLES.translation[f"{wc.sample}_2.fastq.gz"],
    output:
        fq1=temp("capcruncher_output/interim/fastq/{sample}_1.fastq.gz"),
        fq2=temp("capcruncher_output/interim/fastq/{sample}_2.fastq.gz"),
    log:
        "capcruncher_output/logs/fastq_rename/{sample}.log",
    shell:
        """
        ln -s $(realpath {input.fq1}) {output.fq1} &&
        ln -s $(realpath {input.fq2}) {output.fq2}
        """


checkpoint split:
    input:
        fq1=rules.fastq_rename.output.fq1,
        fq2=rules.fastq_rename.output.fq2,
    output:
        directory("capcruncher_output/interim/fastq/split/{sample}"),
    threads: 4
    resources:
        mem_mb=1000,
        time="0-03:00:00",
    params:
        prefix="capcruncher_output/interim/fastq/split/{sample}/{sample}",
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


if not CAPCRUNCHER_TOOLS:

    rule deduplication_parse:
        input:
            fq1="capcruncher_output/interim/fastq/split/{sample}/{sample}_part{part}_1.fastq.gz",
            fq2="capcruncher_output/interim/fastq/split/{sample}/{sample}_part{part}_2.fastq.gz",
        output:
            temp(
                "capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}_{part}.pkl"
            ),
        resources:
            mem_mb=2000,
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
            temp("capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}.pkl"),
        log:
            "capcruncher_output/logs/deduplication_fastq/identify/{sample}.log",
        resources:
            mem_mb=5000,
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
            fq1="capcruncher_output/interim/fastq/split/{sample}/{sample}_part{part}_1.fastq.gz",
            fq2="capcruncher_output/interim/fastq/split/{sample}/{sample}_part{part}_2.fastq.gz",
            ids_duplicated="capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}.pkl",
        output:
            fq1=temp(
                "capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}_part{part}_1.fastq.gz"
            ),
            fq2=temp(
                "capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}_part{part}_2.fastq.gz"
            ),
            stats=temp(
                "capcruncher_output/interim/statistics/deduplication/data/{sample}_part{part}.deduplication.csv"
            ),
        params:
            prefix_fastq="capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}_part{part}",
            prefix_stats="capcruncher_output/interim/statistics/deduplication/data/{sample}_part{part}",
        resources:
            mem_mb=2000,
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
            fq1="capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}_part{part}_1.fastq.gz",
            fq2="capcruncher_output/interim/fastq/deduplicated/{sample}/{sample}_part{part}_2.fastq.gz",
        output:
            trimmed1=temp(
                "capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_1.fastq.gz"
            ),
            trimmed2=temp(
                "capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_2.fastq.gz"
            ),
        params:
            outdir="capcruncher_output/interim/fastq/trimmed/{sample}/",
        threads: 4
        resources:
            mem_mb=2000,
        log:
            "capcruncher_output/logs/trimming/{sample}_{part}.log",
        shell:
            """
            trim_galore --cores {threads} --trim-n --paired --output_dir {params.outdir} {input.fq1} {input.fq2} >> {log} 2>&1 &&
            mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_1_val_1.fq.gz {output.trimmed1} &&
            mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_2_val_2.fq.gz {output.trimmed2}
            """

else:

    checkpoint deduplication:
        input:
            unpack(get_fastq_split_1),
        output:
            fastq_dir=directory(
                "capcruncher_output/interim/fastq/deduplicated/{sample}/"
            ),
            stats="capcruncher_output/interim/statistics/deduplication/data/{sample}.deduplication.csv",
        params:
            prefix_fastq="capcruncher_output/interim/fastq/deduplicated/{sample}/",
            prefix_stats="capcruncher_output/interim/statistics/deduplication/data/{sample}/",
        log:
            "capcruncher_output/logs/deduplication_fastq/{sample}.log",
        threads: workflow.cores * 0.5
        resources:
            mem_mb=lambda wildcards, attempt: 2000 * 2**attempt,
        shell:
            """
            mkdir -p {params.prefix_stats} &&
            capcruncher-tools fastq-deduplicate -1 {input.fq1} -2 {input.fq2} -o {params.prefix_fastq} --statistics {output.stats} --sample-name {wildcards.sample} > {log} 2>&1
            """

    rule trim:
        input:
            unpack(get_deduplicated_fastq_pair),
        output:
            trimmed1=temp(
                "capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_1.fastq.gz"
            ),
            trimmed2=temp(
                "capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_2.fastq.gz"
            ),
        params:
            outdir="capcruncher_output/interim/fastq/trimmed/{sample}/",
        threads: 4
        resources:
            mem_mb=2000,
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
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz"
        ),
        pe1=temp(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_1.fastq.gz"
        ),
        pe2=temp(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_2.fastq.gz"
        ),
        hist=temp(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.hist"
        ),
        histogram=temp(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.histogram"
        ),
    params:
        outdir="capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}",
    threads: 4
    resources:
        mem_mb=1000,
    log:
        "capcruncher_output/logs/flash/{sample}_{part}.log",
    shell:
        """
        flash {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z --compress-prog-args pigz > {log} 2>&1
        """


checkpoint rebalance_partitions_combined:
    input:
        flashed=lambda wc: get_flashed_fastq(wc),
    output:
        directory("capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/"),
        touch(
            "capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/.complete.sentinel"
        ),
    params:
        prefix=lambda wildcards, output: pathlib.Path(output[0]) / wildcards.sample,
        suffix=lambda wc: f"_flashed",
        fq=lambda wc: ",".join(get_flashed_fastq(wc)),
        n_reads=str(config["split"].get("n_reads", 1e6) * 4),
    log:
        "capcruncher_output/logs/rebalance_partitions/{sample}_flashed.log",
    threads: 4
    resources:
        mem_mb=1000,
    shell:
        """
        mkdir -p {output[0]} &&
        capcruncher \
        fastq \
        split \
        {params.fq} \
        -m \
        unix \
        -o \
        {params.prefix} \
        -n \
        {params.n_reads} \
        --gzip \
        -p \
        {threads} \
        --suffix \
        {params.suffix} \
        > {log} 2>&1 &&
        touch {output[1]}
        """


checkpoint rebalance_partitions_pe:
    input:
        fq=get_pe_fastq,
    output:
        directory("capcruncher_output/interim/fastq/rebalanced/{sample}/pe"),
        touch(
            "capcruncher_output/interim/fastq/rebalanced/{sample}/pe/.complete.sentinel"
        ),
    params:
        prefix=lambda wildcards, output: pathlib.Path(output[0]) / wildcards.sample,
        suffix=lambda wc: f"_pe",
        n_reads=str((config["split"].get("n_reads", 1e6) * 4 // 2)),
        fq1=lambda wc: ",".join(separate_pe_fastq(wc)[1]),
        fq2=lambda wc: ",".join(separate_pe_fastq(wc)[2]),
    log:
        "capcruncher_output/logs/rebalance_partitions/{sample}_pe.log",
    threads: 4
    resources:
        mem_mb=1000,
    shell:
        """
        mkdir -p {output[0]} &&
        capcruncher \
        fastq \
        split \
        {params.fq1} \
        {params.fq2} \
        -m \
        unix \
        -o \
        {params.prefix} \
        -n \
        {params.n_reads} \
        --gzip \
        -p \
        {threads} \
        --suffix \
        {params.suffix} \
        > {log} 2>&1 &&
        touch {output[1]}
        """


rule digest_flashed_combined:
    input:
        flashed="capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/{sample}_part{part}_flashed_1.fastq.gz",
        sentinel="capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/.complete.sentinel",
    output:
        digested=temp(
            "capcruncher_output/interim/fastq/digested/{sample}/{sample}_part{part}_flashed.fastq.gz"
        ),
        stats_read=temp(
            "capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_flashed.digestion.read.summary.csv"
        ),
        stats_unfiltered=temp(
            "capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_flashed.digestion.unfiltered.histogram.csv"
        ),
        stats_filtered=temp(
            "capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_flashed.digestion.filtered.histogram.csv"
        ),
    params:
        prefix_stats="capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_flashed",
        restriction_site=config["analysis"]["restriction_enzyme"],
    threads: 4
    resources:
        mem_mb=2000,
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
        pe1="capcruncher_output/interim/fastq/rebalanced/{sample}/pe/{sample}_part{part}_pe_1.fastq.gz",
        pe2="capcruncher_output/interim/fastq/rebalanced/{sample}/pe/{sample}_part{part}_pe_2.fastq.gz",
        sentinel="capcruncher_output/interim/fastq/rebalanced/{sample}/pe/.complete.sentinel",
    output:
        digested=temp(
            "capcruncher_output/interim/fastq/digested/{sample}/{sample}_part{part}_pe.fastq.gz"
        ),
        stats_read=temp(
            "capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_pe.digestion.read.summary.csv"
        ),
        stats_unfiltered=temp(
            "capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_pe.digestion.unfiltered.histogram.csv"
        ),
        stats_filtered=temp(
            "capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_pe.digestion.filtered.histogram.csv"
        ),
    params:
        prefix_stats="capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_pe",
        restriction_site=config["analysis"]["restriction_enzyme"],
    threads: 4
    resources:
        mem_mb=2000,
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
