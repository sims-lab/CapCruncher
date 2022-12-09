import os
import pathlib 
import utils

# def collate_fastq_files_for_splitting(wc):
#     fqs = sorted(list(pathlib.Path(".").glob("*.fastq*")))
#     regex = "(.*?)_[12].fastq.*"
#     fqs_grouped = utils.group_files_by_regex(fqs, regex)

#     for file_group in fqs_grouped.iteritems():

#         if wc.sample in file_group[0]:
#             try:
#                 return {"fastq_1": str(file_group[1][0]), 
#                         "fastq_2": str(file_group[1][1])}
#             except IndexError:
#                 return {}

# def aggregate_split_fastq_names(wc):
#     checkpoint_output = checkpoints.split_fastqs.get(**wc).output[0]
#     sample, part, read, compression = glob_wildcards(os.path.join(checkpoint_output, "{sample}_part{part,\d+}_{read,[12]}.{compression,fastq.gz|fastq}"))
#     file_names = expand(os.path.join(checkpoint_output, "{sample}_part{part}_{read}.{compression}"), sample=sample, part=part, read=read, compression=compression)
#     return file_names

# def collate_fastq_files_for_trimming(wc):

#     fqs = aggregate_split_fastq_names(wc)
#     regex = ".*/(.*?)_part(\d+)_[12].fastq.*"
#     fqs_grouped = utils.group_files_by_regex(fqs, regex)

#     for file_group in fqs_grouped.iteritems():

#         if wc.sample in file_group[0] and wc.part in file_group[0]:
#             try:
#                 return {"fastq_1": str(file_group[1][0]), 
#                         "fastq_2": str(file_group[1][1])}
#             except IndexError:
#                 return {}

# def collate_hashed_reads_for_identification(wc):
#     files = list(pathlib.Path("capcruncher_preprocessing/deduplicated/duplicated_ids/").glob("*.pkl"))
#     return [str(fn) for fn in files if wc.sample in str(fn)]

# def get_paired_fastq_files(wc):

#     return 
    
#     # df = fastq_samples.design.loc[lambda df: df["sample"] == wc.sample]
#     # if not df.empty:
#     #     return {"fq1": df.iloc[0].loc["fq1"], 
#     #             "fq2": df.iloc[0].loc["fq2"]
#     #             }
#     # else:
#     #     return {"fq1": "",
#     #             "fq2": ""}

def get_partitions(wc):
    checkpoint_output = checkpoints.split.get(**wc).output[0]
    parts = glob_wildcards(os.path.join(checkpoint_output, "{sample}_part{part}_{read}.fastq.gz")).part
    return set(parts)

def get_hashed_reads(wc):
    return expand("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_{part}.pkl",
                  sample=wc.sample,
                  part=get_partitions(wc))

def aggregate(wc):
    # files =  expand("capcruncher_preprocessing/flashed/{sample}/{sample}_part{part}_{read}.fastq.gz",
    #               sample=wc.sample,
    #               part=get_partitions(wc),
    #               read=[1,2])
    files =  expand("capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz",
                sample=wc.sample,
                part=get_partitions(wc))
    assert files
    return files



checkpoint split:
    input:
        fq1 = "{sample}_1.fastq.gz",
        fq2 = "{sample}_2.fastq.gz"
    output:
        directory("capcruncher_preprocessing/01_split/{sample}")
    threads:
        2
    params:
        prefix = "capcruncher_preprocessing/01_split/{sample}/{sample}",
        n_reads = str(config["split"].get("n_reads", 1e6)),
        gzip = "--gzip" if COMPRESS_FASTQ else "--no-gzip",
        compression_level = f"--compression_level {config['pipeline'].get('compression', 0)}" if COMPRESS_FASTQ else ""    
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
        fq1 = "capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2 = "capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_{part}.pkl")
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
        hashed_reads = get_hashed_reads
    output:
        temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}.pkl")
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
        fq1 = "capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2 = "capcruncher_preprocessing/01_split/{sample}/{sample}_part{part}_2.fastq.gz",
        ids_duplicated = "capcruncher_preprocessing/02_deduplicated/{sample}/{sample}.pkl",
    output:
        fq1 = temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_1.fastq.gz"),
        fq2 = temp("capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_2.fastq.gz"),
        stats = temp("capcruncher_statistics/deduplication/data/{sample}_part{part}.csv")
    params:
        prefix_fastq = "capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}",
        prefix_stats = "capcruncher_statistics/deduplication/data/{sample}_part{part}",
    threads:
        4
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
    # Trim reads using trimgalore
    input:
        fq1 = "capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2 = "capcruncher_preprocessing/02_deduplicated/{sample}/{sample}_part{part}_2.fastq.gz",
        
    output:
        trimmed1=temp("capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_1.fastq.gz"),
        trimmed2=temp("capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_2.fastq.gz"),
    params:
        outdir = "capcruncher_preprocessing/03_trimmed/{sample}/"
    threads:
        4
    log:
        "logs/trimming/{sample}_{part}.log"
    shell: 
        """
           trim_galore --cores {threads} --trim-n --paired --output_dir {params.outdir} {input.fq1} {input.fq2} >> {log} 2>&1 &&
           mv capcruncher_preprocessing/trimmed/{wildcards.sample}/{wildcards.sample}_part{wildcards.part}_1_val_1.fq.gz {output.trimmed1} &&
           mv capcruncher_preprocessing/trimmed/{wildcards.sample}/{wildcards.sample}_part{wildcards.part}_2_val_2.fq.gz {output.trimmed2}
        """

rule flash:
    input:
        fq1 = "capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2 = "capcruncher_preprocessing/03_trimmed/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        flashed = "capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz",
        pe1 = "capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_1.fastq.gz",
        pe2 = "capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}.notCombined_2.fastq.gz",
    params:
        outdir = "capcruncher_preprocessing/04_flashed/{sample}/{sample}_part{part}"
    threads:
        8
    log:
        "logs/flash/{sample}_{part}.log"
    shell:
        """
        flash {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z --compress-prog-args pigz > {log} 2>&1
        """


rule done:
    input:
        unpack(aggregate)
    output:
        touch("{sample}.done")




# rule deduplication_parse:
#     input:
#         fq1 = "capcruncher_preprocessing/trimmed/{sample}/{sample}_part{part}_1.fastq.gz",
#         fq2 = "capcruncher_preprocessing/trimmed/{sample}/{sample}_part{part}_2.fastq.gz",
#     output:
#         temp("capcruncher_preprocessing/deduplicated/duplicated_ids/{sample}_{part}.pkl")
#     shell:
#         """
#         capcruncher \
#         fastq \
#         deduplicate \
#         parse \
#         {input.fq1} \
#         {input.fq2} \
#         -o \
#         {output}
#         """

# rule deduplication_identify:
#     input:
#         hashed_reads = expand("capcruncher_preprocessing/deduplicated/duplicated_ids/{sample}/{sample}_{part}.pkl"),
#     output:
#         temp("capcruncher_preprocessing/deduplicated/duplicated_ids/{sample}/{sample}.pkl")
#     shell:
#         """
#         capcruncher \
#         fastq \
#         deduplicate \
#         identify \
#         {input.hashed_reads} \
#         -o \
#         {output}
#         """

# # rule test_stop:
# #     input:
        
# #     output:
# #         touch("{sample}.done")
    




