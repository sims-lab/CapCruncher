import capcruncher.pipeline.utils
from typing import Literal


def get_rebalanced_parts(
    wildcards, combined: Literal["flashed", "pe"] = None, **kwargs
):
    combined = combined or wildcards.combined
    import pathlib
    import re

    parts = dict()
    outdirs = dict(
        flashed=checkpoints.rebalance_partitions_combined.get(
            **{**wildcards, **kwargs}
        ).output[0],
        pe=checkpoints.rebalance_partitions_pe.get(**{**wildcards, **kwargs}).output[0],
    )

    for combined_type in ["flashed", "pe"]:
        fq_files = pathlib.Path(outdirs[combined_type]).glob("*.fastq.gz")
        parts[combined_type] = list(
            sorted(
                set(
                    [
                        int(re.search(r"part(\d+)", f.name).group(1))
                        for f in fq_files
                        if re.search(r"part(\d+)", f.name)
                    ]
                )
            )
        )

    if combined == "flashed":
        return parts["flashed"]
    else:
        return parts["pe"]


def get_rebalanced_bam(wildcards):
    bam = []
    for combined_type in ["flashed", "pe"]:
        for part in get_rebalanced_parts(wildcards, combined_type):
            bam.append(
                f"capcruncher_output/interim/aligned/{wildcards.sample}/{wildcards.sample}_part{part}_{combined_type}.sorted.bam"
            )

    return bam


rule align_bowtie2:
    input:
        fastq="capcruncher_output/interim/fastq/digested/{sample}/{sample}_part{part}_{combined}.fastq.gz",
    output:
        bam="capcruncher_output/interim/aligned/{sample}/{sample}_part{part}_{combined,(flashed|pe)}.bam",
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
        bam="capcruncher_output/interim/aligned/{sample}/{sample}_part{part}_{combined}.sorted.bam",
    threads: 4
    log:
        "capcruncher_output/logs/align/{sample}_{part}_{combined}_sort.log",
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.bam} 2> {log}
        """


rule merge_bam_partitions:
    input:
        bam=get_rebalanced_bam,
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
