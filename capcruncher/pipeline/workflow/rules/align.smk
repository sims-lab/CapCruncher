rule align_bowtie2:
    input:
        fastq="capcruncher_output/fastq/digested/{sample}/{sample}_part{part}_{combined}.fastq.gz",
    output:
        bam=temp(
            "capcruncher_output/aligned/{sample}/{sample}_part{part}_{combined,(flashed|pe)}.bam"
        ),
    resources:
        mem_mb=4000,
    params:
        aligner=config["align"]["aligner"],
        index_flag=config["align"].get("index_flag", ""),
        indices=config["genome"]["aligner_index"],
        options=config["align"].get("options", ""),
    threads: 4
    log:
        "capcruncher_output/logs/align/{sample}_{part}_{combined}.log",
    shell:
        """
        {params.aligner} {params.index_flag} {params.indices} {params.options} -p {threads} {input.fastq} 2> {log} |
        samtools view -bS - > {output.bam}
        """


rule sort_bam_partitions:
    input:
        bam=rules.align_bowtie2.output.bam,
    output:
        bam=temp(
            "capcruncher_output/aligned/{sample}/{sample}_part{part}_{combined}.sorted.bam"
        ),
    threads: 4
    log:
        "capcruncher_output/logs/align/{sample}_{part}_{combined}_sort.log",
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam} 2> {log}
        """


rule merge_bam_partitions:
    input:
        bam=lambda wc: expand(
            "capcruncher_output/aligned/{{sample}}/{{sample}}_part{part}_{combined}.sorted.bam",
            part=get_parts(wc),
            combined=["flashed", "pe"],
        ),
    output:
        bam="capcruncher_output/aligned/{sample}.bam",
    shell:
        """
        samtools merge -o {output.bam} {input.bam}
        """


rule index_bam:
    input:
        bam=rules.merge_bam_partitions.output.bam,
    output:
        bam="capcruncher_output/aligned/{sample}.bam.bai",
    shell:
        """
        samtools index {input.bam}
        """
