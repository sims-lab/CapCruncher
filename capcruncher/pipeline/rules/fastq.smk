import os
import pathlib
import utils


def get_fastq_partition_numbers_for_sample(wc):

    checkpoint_output = checkpoints.split.get(**wc).output[0]
    sample_name = wc.sample
    parts = glob_wildcards(
        os.path.join(
            checkpoint_output, "".join([sample_name, "_part{part}_{read}.fastq.gz"])
        )
    ).part

    return set(parts)


def get_hashed_reads(wc):
    return expand(
        "capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_{part}.pkl",
        sample=wc.sample,
        part=get_fastq_partition_numbers_for_sample(wc),
    )


# def aggregate(wc):
#     files = expand(
#         "capcruncher_preprocessing/05_digested/{sample}/{sample}_part{part}_{combined}.fastq.gz",
#         sample=wc.sample,
#         part=get_fastq_partition_numbers_for_sample(wc),
#         combined=["flashed", "pe"],
#     )
#     assert files
#     return files


checkpoint split:
    input:
        fq1="{sample}_1.fastq.gz",
        fq2="{sample}_2.fastq.gz",
    output:
        directory("capcruncher_preprocessing/01_split/{sample}"),
    threads: 2
    params:
        prefix="capcruncher_preprocessing/01_split/{sample}/{sample}",
        n_reads=str(config["split"].get("n_reads", 1e6)),
        gzip="--gzip" if COMPRESS_FASTQ else "--no-gzip",
        compression_level=f"--compression_level {config['pipeline'].get('compression', 0)}"
        if COMPRESS_FASTQ
        else "",
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
        {threads}
        """


rule deduplication_parse:
    input:
        fq1="capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_{part}.pkl"),
    shell:
        """
        capcruncher \
        fastq \
        deduplicate \
        parse \
        {input.fq1} \
        {input.fq2} \
        -o \
        {output}
        """


rule deduplication_identify:
    input:
        hashed_reads=get_hashed_reads,
    output:
        temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}.pkl"),
    shell:
        """
        capcruncher \
        fastq \
        deduplicate \
        identify \
        {input.hashed_reads} \
        -o \
        {output}
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
        stats=temp(
            "capcruncher_statistics/deduplication/data/{sample}_part{part}.deduplication.csv"
        ),
    params:
        prefix_fastq="capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}",
        prefix_stats="capcruncher_statistics/deduplication/data/{sample}_part{part}",
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
        --gzip
        """


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
        flashed="capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz",
        pe1="capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_1.fastq.gz",
        pe2="capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_2.fastq.gz",
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
        digested="capcruncher_preprocessing/05_digested/{sample}/{sample}_part{part}_flashed.fastq.gz",
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
        digested="capcruncher_preprocessing/05_digested/{sample}/{sample}_part{part}_pe.fastq.gz",
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


# rule done:
#     input:
#         unpack(aggregate),
#     output:
#         touch("{sample}.done"),
