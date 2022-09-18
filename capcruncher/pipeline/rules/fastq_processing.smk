import os
import pathlib 
import utils

def collate_fastq_files_for_splitting(wc):
    fqs = sorted(list(pathlib.Path(".").glob("*.fastq*")))
    regex = "(.*?)_[12].fastq.*"
    fqs_grouped = utils.group_files_by_regex(fqs, regex)

    for file_group in fqs_grouped.iteritems():

        if wc.sample in file_group[0]:
            try:
                return {"fastq_1": str(file_group[1][0]), 
                        "fastq_2": str(file_group[1][1])}
            except IndexError:
                return {}

def aggregate_split_fastq_names(wc):
    checkpoint_output = checkpoints.split_fastqs.get(**wc).output[0]
    sample, part, read, compression = glob_wildcards(os.path.join(checkpoint_output, "{sample}_part{part,\d+}_{read,[12]}.{compression,fastq.gz|fastq}"))
    file_names = expand(os.path.join(checkpoint_output, "{sample}_part{part}_{read}.{compression}"), sample=sample, part=part, read=read, compression=compression)
    return file_names

def collate_fastq_files_for_trimming(wc):

    fqs = aggregate_split_fastq_names(wc)
    regex = ".*/(.*?)_part(\d+)_[12].fastq.*"
    fqs_grouped = utils.group_files_by_regex(fqs, regex)

    for file_group in fqs_grouped.iteritems():

        if wc.sample in file_group[0] and wc.part in file_group[0]:
            try:
                return {"fastq_1": str(file_group[1][0]), 
                        "fastq_2": str(file_group[1][1])}
            except IndexError:
                return {}

def collate_hashed_reads_for_identification(wc):
    files = list(pathlib.Path("capcruncher_preprocessing/deduplicated/duplicated_ids/").glob("*.pkl"))
    return [str(fn) for fn in files if wc.sample in str(fn)]



checkpoint split_fastqs:
    input:
        unpack(collate_fastq_files_for_splitting)
    output:
        directory("capcruncher_preprocessing/split/{sample}")
    threads:
        2
    params:
        n_reads = str(config["split"].get("n_reads", 1e6)),
        gzip = "--gzip" if COMPRESS_FASTQ else "--no-gzip",
        compression_level = f"--compression_level {config['pipeline'].get('compression', 0)}" if COMPRESS_FASTQ else ""    
    shell:
        """
        mkdir {output} && \
        capcruncher \
        fastq \
        split \
        {input.fastq_1} \
        {input.fastq_2} \
        -m \
        unix \
        -o \
        capcruncher_preprocessing/split/{wildcards.sample}/{wildcards.sample} \
        -n \
        {params.n_reads} \
        {params.gzip} \
        {params.compression_level} \
        -p \
        {threads}
        """

rule trim:
    # Trim reads using trimgalore
    input:
        unpack(collate_fastq_files_for_trimming),
    output:
        trimmed1=temp("capcruncher_preprocessing/trimmed/{sample}/{sample}_part{part}_1.fastq.gz"),
        trimmed2=temp("capcruncher_preprocessing/trimmed/{sample}/{sample}_part{part}_2.fastq.gz"),
    threads:
        4
    log:
        "logs/trimming/{sample}_{part}.log"
    shell: 
        """
           trim_galore --cores {threads} --trim-n --paired --output_dir trimmed {input.fastq_1} {input.fastq_2} >> {log} 2>&1 &&
           mv capcruncher_preprocessing/trimmed/{wildcards.sample}/{wildcards.sample}_part{part}_1_val_1.fq.gz {output.trimmed1} &&
           mv capcruncher_preprocessing/trimmed/{wildcards.sample}/{wildcards.sample}_part{part}_2_val_2.fq.gz {output.trimmed2}
        """

rule deduplication_parse:
    input:
        fq1 = "capcruncher_preprocessing/trimmed/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2 = "capcruncher_preprocessing/trimmed/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        temp("capcruncher_preprocessing/deduplicated/duplicated_ids/{sample}_{part}.pkl")
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
        hashed_reads = collate_hashed_reads_for_identification
    output:
        temp("capcruncher_preprocessing/deduplicated/duplicated_ids/{sample}/{sample}.pkl")
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

rule test_stop:
    input: 
        "capcruncher_preprocessing/deduplicated/duplicated_ids/{sample}/{sample}.pkl"
    output:
        touch("{sample}.done")
    




