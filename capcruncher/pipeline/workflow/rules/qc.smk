import os
import pathlib
import capcruncher.pipeline.utils


rule fastqc:
    input:
        "capcruncher_output/interim/fastq/{sample}_{read}.fastq.gz",
    output:
        html="capcruncher_output/interim/qc/fastqc/{sample}_{read}_fastqc.html",
        zip="capcruncher_output/interim/qc/fastqc/{sample}_{read}_fastqc.zip",
    log:
        "capcruncher_output/logs/fastqc/{sample}_{read}.log",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: scale_memory(1, attempt),
    params:
        extra="--quiet",
        outdir=lambda wc, output: pathlib.Path(output.html).parent,
        memory=1024,
    shell:
        """
        mkdir -p {params.outdir} \
            && fastqc \
                --threads {threads} \
                --memory {params.memory} \
                {params.extra} \
                --outdir {params.outdir} \
                {input} \
                >{log} 2>&1
        """


rule samtools_stats:
    input:
        bam="capcruncher_output/results/{sample}/{sample}.bam",
        bai="capcruncher_output/results/{sample}/{sample}.bam.bai",
    output:
        stats=temp("capcruncher_output/interim/qc/alignment_raw/{sample}.txt"),
    log:
        "capcruncher_output/logs/samtools_stats/{sample}.log",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: scale_memory(1, attempt),
    shell:
        """
        samtools stats {input.bam} >{output.stats} 2>{log}
        """


rule multiqc_report:
    input:
        expand(
            "capcruncher_output/interim/qc/fastqc/{sample}_{read}_fastqc.html",
            sample=SAMPLE_NAMES,
            read=[1, 2],
        ),
        expand(
            "capcruncher_output/interim/qc/alignment_raw/{sample}.txt",
            sample=SAMPLE_NAMES,
        ),
    output:
        "capcruncher_output/results/full_qc_report.html",
    log:
        "capcruncher_output/logs/multiqc.log",
    resources:
        mem=lambda wildcards, attempt: scale_memory(1, attempt),
    params:
        outdir=lambda wc, output: str(pathlib.Path(output[0]).parent),
        dir_analysis=lambda wc, input: str(pathlib.Path(input[0]).parents[1]),
    shell:
        "multiqc -o {params.outdir} {params.dir_analysis} -n full_qc_report.html --force > {log} 2>&1"


rule multiqc_full:
    input:
        report="capcruncher_output/results/full_qc_report.html",
    output:
        report="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_report.html",
        trimming_data="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_cutadapt.txt",
        flash_data="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_flash_combo_stats.txt",
        bowtie2_data="capcruncher_output/interim/statistics/multiqc_full_data/multiqc_data/multiqc_bowtie2.txt",
    log:
        "capcruncher_output/logs/multiqc_full.log",
    shell:
        "multiqc capcruncher_output/ --outdir capcruncher_output/interim/statistics/multiqc_full_data --force > {log} 2>&1"
