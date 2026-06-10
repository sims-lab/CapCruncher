rule align_bowtie2:
    input:
        fastq="capcruncher_output/interim/fastq/digested/{sample}/{sample}_part{part}_{combined}.fastq.gz",
    output:
        bam=temp(
            "capcruncher_output/interim/aligned/{sample}/{sample}_part{part}_{combined,(flashed|pe)}.bam"
        ),
    log:
        "capcruncher_output/logs/align/{sample}_{part}_{combined}.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(4, attempt),
    params:
        aligner=config["align"]["aligner"],
        index_flag=config["align"].get("index_flag", ""),
        indices=config["genome"]["aligner_index"],
        options=config["align"].get("options", ""),
    shell:
        """
        {params.aligner} {params.index_flag} {params.indices} {params.options} -p {threads} {input.fastq} 2>{log} \
            | samtools view -bS - >{output.bam}
        """


rule sort_bam_partitions:
    input:
        bam=rules.align_bowtie2.output.bam,
    output:
        bam=temp(
            "capcruncher_output/interim/aligned/{sample}/{sample}_part{part}_{combined}.sorted.bam"
        ),
    log:
        "capcruncher_output/logs/align/{sample}_{part}_{combined}_sort.log",
    threads: 4
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam} 2>{log}
        """


rule merge_bam_partitions:
    input:
        bam=get_rebalanced_bam,
    output:
        bam="capcruncher_output/results/{sample}/{sample}.bam",
    log:
        "capcruncher_output/logs/merge_bam_partitions/{sample}.log",
    shell:
        """
        samtools merge {output.bam} {input.bam} >{log} 2>&1
        """


rule index_bam:
    input:
        bam="capcruncher_output/results/{sample}/{sample}.bam",
    output:
        bam="capcruncher_output/results/{sample}/{sample}.bam.bai",
    log:
        "capcruncher_output/logs/index_bam/{sample}.log",
    shell:
        """
        samtools index {input.bam} >{log} 2>&1
        """
