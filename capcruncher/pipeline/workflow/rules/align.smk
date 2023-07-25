import capcruncher.pipeline.utils


def get_rebalanced_parts(wildcards, combined: Literal["flashed", "pe"] = None):

    # Confirm that the checkpoint has been run
    checkpoints.rebalance_partitions_combined.get(**wildcards)

    combined = combined or wildcards.combined
    with open("capcruncher_output/resources/rebalanced/{wildcards.sample}.json") as f:
        rebalanced = json.load(f)

    if wildcards.combined == "flashed":
        return rebalanced["flashed"]
    else:
        return rebalanced["pe"]


def get_rebalanced_bam(wildcards):
    bam = []
    for combined_type in ["flashed", "pe"]:
        for part in get_rebalanced_parts(wildcards, combined_type):
            bam.append(
                f"capcruncher_output/interim/aligned/{wildcards.sample}/{wildcards.sample}_part{part}_{combined_type}.sorted.bam"
            )


rule align_bowtie2:
    input:
        fastq="capcruncher_output/interim/fastq/digested/{sample}/{sample}_part{part}_{combined}.fastq.gz",
    output:
        bam=temp(
            "capcruncher_output/interim/aligned/{sample}/{sample}_part{part}_{combined,(flashed|pe)}.bam"
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
            "capcruncher_output/interim/aligned/{sample}/{sample}_part{part}_{combined}.sorted.bam"
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
        bam=capcruncher.pipeline.utils.get_rebalanced_bam,
        n_parts="capcruncher_output/resources/rebalanced/{sample}.json",
    output:
        bam="capcruncher_output/results/{sample}/{sample}.bam",
    shell:
        """
        samtools merge -o {output.bam} {input.bam}
        """


rule index_bam:
    input:
        bam="capcruncher_output/results/{sample}/{sample}.bam",
    output:
        bam="capcruncher_output/results/{sample}/{sample}.bam.bai",
    shell:
        """
        samtools index {input.bam}
        """
