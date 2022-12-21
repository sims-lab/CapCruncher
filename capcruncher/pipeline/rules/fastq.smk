import os
import pathlib
import utils
import json


def get_fastq_partition_numbers_for_sample(wc, sample_name=None):

    if not sample_name:
        sample_name = wc.sample

    n_partitions_path = f"capcruncher_resources/n_partitions/{sample_name}.json"

    if os.path.exists(n_partitions_path):
        with open(n_partitions_path, "r") as f:
            n_partitions = json.load(f)
    else:
        checkpoint_output = checkpoints.split.get(sample=sample_name).output[0]
        n_partitions = glob_wildcards(
            os.path.join(
                checkpoint_output, "".join([sample_name, "_part{part}_{read}.fastq.gz"])
            )
        ).part

        n_partitions = list(set(n_partitions))

        if not os.path.exists("capcruncher_resources/n_partitions"):
            os.makedirs("capcruncher_resources/n_partitions")

        with open(n_partitions_path, "w") as f:
            json.dump(n_partitions, f)

    return n_partitions


checkpoint split:
    input:
        fq1="{sample}_1.fastq.gz",
        fq2="{sample}_2.fastq.gz",
    output:
        directory(temp("capcruncher_preprocessing/01_split/{sample}")),
    threads: 2
    params:
        prefix="capcruncher_preprocessing/01_split/{sample}/{sample}",
        n_reads=str(config["split"].get("n_reads", 1e6)),
        gzip="--gzip" if COMPRESS_FASTQ else "--no-gzip",
        compression_level=f"--compression_level {config['pipeline'].get('compression', 0)}"
        if COMPRESS_FASTQ
        else "",
    log:
        "logs/split/{sample}.log",
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
        {params.gzip} \
        {params.compression_level} \
        -p \
        {threads} \
        > {log} 2>&1
        """


rule deduplication_parse:
    input:
        fq1="capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_{part}.pkl"),
    log:
        "logs/deduplication_fastq/parse/{sample}_part{part}.log",
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
        hashed_reads=lambda wc: expand(
            "capcruncher_preprocessing/02_deduplicated/{{sample}}/{{sample}}_{part}.pkl",
            part=get_fastq_partition_numbers_for_sample(wc),
        ),
    output:
        temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}.pkl"),
    log:
        "logs/deduplication_fastq/identify/{sample}.log",
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
        fq1="capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_2.fastq.gz",
        ids_duplicated="capcruncher_preprocessing/02_deduplicated/{sample}/{sample}.pkl",
    output:
        fq1=temp(
            "capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_1.fastq.gz"
        ),
        fq2=temp(
            "capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_2.fastq.gz"
        ),
        stats="capcruncher_statistics/deduplication/data/{sample}_part{part}.deduplication.csv",
    params:
        prefix_fastq="capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}",
        prefix_stats="capcruncher_statistics/deduplication/data/{sample}_part{part}",
    log:
        "logs/deduplication_fastq/remove/{sample}_part{part}.log",
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


# rule deduplication_finished:
#     input:
#         deduplication_performed=lambda wc: expand(
#             "capcruncher_statistics/deduplication/data/{{sample}}_part{part}.deduplication.csv",
#             part=get_fastq_partition_numbers_for_sample(wc),
#         ),
#     output:
#         "capcruncher_reseources/sentinels/{sample}.deduplication_finished.sentinel",
#     shell:
#         """
#         touch {output}
#         """


# rule remove_split_files:
#     input:
#         split_dir="capcruncher_preprocessing/01_split/{sample}/",
#         deduplication_finished=rules.deduplication_finished.output,
#     output:
#         "capcruncher_resources/n_partitions/{sample}_removed.sentinel",
#     log:
#         "logs/remove_split_files/{sample}.log",
#     shell:
#         """
#         rm -f {input.split_dir}/*.fastq.gz && \
#         touch {output} && \
#         echo "Removed split files for {wildcards.sample}" > {log}
#         """


rule trim:
    input:
        fq1="capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        trimmed1=temp(
            "capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_1.fastq.gz"
        ),
        trimmed2=temp(
            "capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_2.fastq.gz"
        ),
    params:
        outdir="capcruncher_preprocessing/03_trimmed/{sample}/",
    threads: 4
    log:
        "logs/trimming/{sample}_{part}.log",
    shell:
        """
           trim_galore --cores {threads} --trim-n --paired --output_dir {params.outdir} {input.fq1} {input.fq2} >> {log} 2>&1 &&
           mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_1_val_1.fq.gz {output.trimmed1} &&
           mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_2_val_2.fq.gz {output.trimmed2}
        """


rule flash:
    input:
        fq1="capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        flashed=temp(
            "capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz"
        ),
        pe1=temp(
            "capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_1.fastq.gz"
        ),
        pe2=temp(
            "capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_2.fastq.gz"
        ),
    params:
        outdir="capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}",
    threads: 8
    log:
        "logs/flash/{sample}_{part}.log",
    shell:
        """
        flash {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z --compress-prog-args pigz > {log} 2>&1
        """


rule digest_flashed_combined:
    input:
        flashed="capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz",
    output:
        digested=temp(
            "capcruncher_preprocessing/05_digested/{sample}/{sample}_part{part}_flashed.fastq.gz"
        ),
        stats_read="capcruncher_statistics/digestion/data/{sample}_part{part}_flashed.digestion.read.summary.csv",
        stats_unfiltered="capcruncher_statistics/digestion/data/{sample}_part{part}_flashed.digestion.unfiltered.histogram.csv",
        stats_filtered="capcruncher_statistics/digestion/data/{sample}_part{part}_flashed.digestion.filtered.histogram.csv",
    params:
        prefix_stats="capcruncher_statistics/digestion/data/{sample}_part{part}_flashed",
        restriction_site=config["analysis"]["restriction_enzyme"],
    threads: 4
    log:
        "logs/digestion/{sample}_{part}.log",
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
        pe1="capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_1.fastq.gz",
        pe2="capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_2.fastq.gz",
    output:
        digested=temp(
            "capcruncher_preprocessing/05_digested/{sample}/{sample}_part{part}_pe.fastq.gz"
        ),
        stats_read="capcruncher_statistics/digestion/data/{sample}_part{part}_pe.digestion.read.summary.csv",
        stats_unfiltered="capcruncher_statistics/digestion/data/{sample}_part{part}_pe.digestion.unfiltered.histogram.csv",
        stats_filtered="capcruncher_statistics/digestion/data/{sample}_part{part}_pe.digestion.filtered.histogram.csv",
    params:
        prefix_stats="capcruncher_statistics/digestion/data/{sample}_part{part}_pe",
        restriction_site=config["analysis"]["restriction_enzyme"],
    threads: 4
    log:
        "logs/digestion/{sample}_{part}.log",
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
